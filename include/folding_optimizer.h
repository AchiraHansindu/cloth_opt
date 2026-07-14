#pragma once

// folding_optimizer.h
//
// Optimization-based cloth folding trajectories for the ClothOpt mass-spring
// engine. Implements FOLDING_SOLUTION_v2: a fold is modeled as a rigid
// rotation of a moving half M about a crease line, executed as per-vertex
// semicircular arcs with cosine easing, driven through the existing
// ClothController position controls.
//
// Geometry (FOLDING_SOLUTION_v2 §2):
//   The sheet lies flat in the xz-plane at working height h (y up). A fold is
//   (crease line, moving set M, stationary set S). The crease is a point a
//   (at y = h) plus a unit in-plane direction dHat; nHat is the unit in-plane
//   normal pointing toward the side that moves. Per vertex:
//       r_v = (p0_v - a) . nHat     (signed distance; r > tol => moving)
//       u_v = (p0_v - a) . dHat
//   M = { r_v > tol }; S = everything else (including crease vertices with
//   |r| <= tol); all of S is pinned (pinning substitutes for table friction
//   in this engine, fact E9). tol = 1e-9.
//
// Trajectory (§5) and timing (§6): see FoldingOptimizer::arc.
//
// Engine facts this module is written against (audited in REPO_REFERENCE.md):
//   E5  position control is pure P: F = gain*(target - pos), clamped to
//       maxForce, applied as velocity += F*dt at assumed unit mass.
//   E6  setTrajectory silently resets gain/maxForce to 1000/100 and returns
//       the LAST waypoint for query times before times[0]; both traps are
//       neutralized in applyFoldingTrajectory (the single choke point).
//   E10 the controlled set and the pinned set must be disjoint (drivers
//       assert this): controlling a pinned vertex accumulates velocity
//       forever while the position never moves.

#include <Eigen/Dense>
#include <cstddef>
#include <vector>

namespace ClothOpt {

class ClothMesh;       // forward declarations; full types needed only in .cpp
class ClothController;

enum class FoldDirection { AlongX, AlongZ };   // crease parallel to this axis; far half moves
enum class Corner { NearLeft, NearRight, FarLeft, FarRight };
// (i,j): NearLeft=(0,0) NearRight=(0,W-1) FarLeft=(H-1,0) FarRight=(H-1,W-1)
enum class GripMode { All, FarEdge, TwoCorners, OneCorner };

struct TrajectoryPoint { size_t vertexIndex; Eigen::Vector3d position; double time; };

struct FoldSpec {
  std::vector<size_t> moving;                  // M
  std::vector<size_t> pinnedSet;               // S incl. crease vertices
  std::vector<Eigen::Vector3d> mirrorTargets;  // ideal (delta = 0) mirror + eps*yHat, same order as moving
  Eigen::Vector3d creasePoint, creaseDir, creaseNormal;   // a, dHat, nHat (nHat -> moving side)
};

struct FoldParams {
  double T = 5.0, gain = 2000.0, maxForce = 400.0, eps = 0.02, delta = 0.0;
  int K = 50;
  GripMode grip = GripMode::All;
  double stiffness = 800.0, bending = 20.0, damping = 0.9, dt = 0.01, settle = 1.5, h = 0.25;
  bool realGravity = false, layerClamp = false;

  // FOLDING_FINAL_FORMULATION §elliptical-family: alpha is the elliptical
  // height ratio scaling ONLY the vertical (yHat) semi-axis of each vertex
  // arc. alpha = 1 is the exact rigid semicircle (zero stretch); alpha < 1
  // lowers the path so the material column (length r) has slack (1-alpha)*r
  // and buckles (less lift, some stretch); alpha > 1 stretches the sheet and
  // is Pareto-dominated. Search allows up to 1.2 so the optimizer confirms
  // rather than assumes the alpha* ~= 1 optimum.
  double alpha = 1.0;        // elliptical height ratio: 1 = rigid arc; <1 lowers path (slack); >1 stretches
  int    substeps = 1;       // integrator substeps per control step (stability for stiff / restored gravity)
  double gravityScale = 1.0; // multiplies restored gravity; continuous replacement for the realGravity boolean
  // gravityScale resolution (backward compat with realGravity): if
  // realGravity is true and gravityScale == 1.0, use gravityScale = N*N (the
  // old E2b behavior); otherwise gravityScale is used directly. See
  // resolveGravityScale() in the driver.
};

struct FoldMetrics {
  double foldErrMean = 0, foldErrMax = 0;      // M1, fraction of L
  double interiorStrainPeak = 0;               // M2, over non-crossing springs, whole rollout
  double touchStrainPeak = 0;                  // peak over interior springs with exactly one endpoint |r| <= tol (crease-touching)
  double crossingMinLenRatio = 1;              // min crossing length / rest (companion of M3)
  int    interpenetrations = 0;                // M6
  double residualKE = 0;                       // M8, mean KE of M over last 1 s of settle
  bool   finite = true;                        // NaN guard
  double gripperEnergy = 0;  // E: sum over steps and gripped vertices of max(0, F . dx), non-regenerative
};

class FoldingOptimizer {
public:
  FoldingOptimizer(int W, int H, double s);

  // Fold geometry per §2. Symmetric fold over Z (crease parallel to X,
  // FoldDirection::AlongX): z_c = (H-1)s/2; a = (0, h, z_c); dHat = (1,0,0);
  // nHat = (0,0,1); mirror (x, h + eps, 2z_c - z). AlongZ swaps roles
  // (x_c = (W-1)s/2, nHat = (1,0,0), mirror x -> 2x_c - x). Convention: the
  // half with the LARGER coordinate moves.
  //
  // Note on the trailing `eps` argument: mirrorTargets are defined (§2) as
  // the ideal delta = 0 mirror PLUS eps*yHat, so the spec must know eps; it
  // is accepted as a defaulted trailing argument (default matches
  // FoldParams::eps) to keep scoring targets consistent with the arcs.
  // `delta` shapes only the arc (r_eff = r - delta, §5) and NEVER the
  // scoring targets; it is accepted here for interface completeness.
  FoldSpec computeFoldSpec(double h, FoldDirection dir, double delta = 0.0,
                           double eps = 0.02) const;                    // symmetric

  // Diagonal fold (§2, general form, supports both diagonals): the crease is
  // the perpendicular bisector of the startCorner-endCorner segment, which
  // for a square is the other diagonal. a = cloth center; nHat =
  // normalize(pos(startCorner) - a) with y zeroed; dHat =
  // normalize(yHat x nHat). Mirror: p' = p - 2 r_p nHat + eps*yHat.
  // Requires W == H (asserted; prints a clear message otherwise).
  FoldSpec computeFoldSpec(double h, Corner startCorner, Corner endCorner,
                           double delta = 0.0, double eps = 0.02) const; // diagonal

  // Which moving vertices the "gripper" drives. Non-gripped moving vertices
  // are left free (neither pinned nor controlled); that is the point of the
  // grip-mode ablation study.
  std::vector<size_t> gripSet(const FoldSpec&, GripMode) const;

  // Semicircular arc about the crease with cosine easing (§5, §6).
  std::vector<Eigen::Vector3d> arc(const Eigen::Vector3d& p0, const FoldSpec&, const FoldParams&,
                                   std::vector<double>& timesOut) const;

  // README-compatible surface: compute spec, take the grip set, generate
  // arcs, return one flat vector sorted by (vertexIndex, time). These do NOT
  // touch the controller.
  std::vector<TrajectoryPoint> optimizeSymmetricFold(const ClothMesh&, FoldDirection, double foldDuration = 5.0,
                                                     FoldParams params = {});
  std::vector<TrajectoryPoint> optimizeDiagonalFold(const ClothMesh&, Corner startCorner, Corner endCorner,
                                                    double foldDuration = 6.0, FoldParams params = {});

  // The single place in the codebase that calls setTrajectory (E6 choke
  // point): re-arms gain/maxForce after every setTrajectory call.
  void applyFoldingTrajectory(ClothController&, const std::vector<TrajectoryPoint>&,
                              double gain, double maxForce);

  // g-fold trajectory (van den Berg / Miller et al. eq. (1)): the gravity-
  // based real-world folding path. A gripped point at fold-line distance y_v
  // follows x = x_v, y = y_b, z = y_v - |y_b| as the baseline y_b sweeps
  // y_v -> -y_v: a 45-degree rise to peak height r above the crease then a
  // mirrored descent, with the ungripped cloth left to hang (needs
  // gravityScale > 1 and substeps to read visibly). Closed-form path — the
  // papers' "optimization" is perception-side, not the trajectory. Drives
  // the OneCorner diagonal grip. See FOLDING_FINAL_FORMULATION §g-fold.
  std::vector<TrajectoryPoint> gFoldTrajectory(const ClothMesh&, Corner startCorner,
                                               Corner endCorner, FoldParams params);

private:
  int W_, H_; double s_;
};

} // namespace ClothOpt

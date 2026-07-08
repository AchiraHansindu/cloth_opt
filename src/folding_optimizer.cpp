// folding_optimizer.cpp — see include/folding_optimizer.h for the overview.
//
// All formulas cite FOLDING_SOLUTION_v2 section numbers. The code is meant to
// be readable as study notes: each core function opens with its derivation.

#include "folding_optimizer.h"
#include "cloth.h"
#include "controller.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <map>

namespace ClothOpt {

namespace {

constexpr double kTol = 1e-9;                       // §2: crease classification tolerance
constexpr double kPi  = 3.14159265358979323846;
const Eigen::Vector3d kYHat(0.0, 1.0, 0.0);

// §2: classify every lattice vertex against the crease frame stored in spec,
// and record the ideal mirror target for each moving vertex.
//   r_v = (p0_v - a) . nHat ;  M = { r_v > tol } ;  S = the rest (incl. crease)
//   mirror(p0) = p0 - 2 r_v nHat + eps*yHat
// The +eps*yHat lift is the target stacking gap: with no self-collision in
// the engine (E8) nothing else keeps the folded flap above the lower half.
// delta NEVER enters here — mirror targets are ALWAYS the ideal delta = 0
// mirror (§2); delta shapes only the arc (§5).
void classifyAndMirror(int W, int H, double s, double h, double eps, FoldSpec& spec) {
  for (int i = 0; i < H; ++i) {
    for (int j = 0; j < W; ++j) {
      const size_t v = static_cast<size_t>(i) * static_cast<size_t>(W) + static_cast<size_t>(j);
      const Eigen::Vector3d p0(j * s, h, i * s);   // createGrid layout: (j*s, y, i*s)
      const double r = (p0 - spec.creasePoint).dot(spec.creaseNormal);
      if (r > kTol) {
        spec.moving.push_back(v);
        spec.mirrorTargets.push_back(p0 - 2.0 * r * spec.creaseNormal + eps * kYHat);
      } else {
        spec.pinnedSet.push_back(v);
      }
    }
  }
}

Eigen::Vector3d cornerPosition(Corner c, int W, int H, double s, double h) {
  switch (c) {
    case Corner::NearLeft:  return Eigen::Vector3d(0.0,         h, 0.0);
    case Corner::NearRight: return Eigen::Vector3d((W - 1) * s, h, 0.0);
    case Corner::FarLeft:   return Eigen::Vector3d(0.0,         h, (H - 1) * s);
    case Corner::FarRight:  return Eigen::Vector3d((W - 1) * s, h, (H - 1) * s);
  }
  return Eigen::Vector3d::Zero();
}

bool cornersAreOpposite(Corner a, Corner b) {
  return (a == Corner::NearLeft  && b == Corner::FarRight) ||
         (a == Corner::FarRight  && b == Corner::NearLeft) ||
         (a == Corner::NearRight && b == Corner::FarLeft)  ||
         (a == Corner::FarLeft   && b == Corner::NearRight);
}

} // namespace

FoldingOptimizer::FoldingOptimizer(int W, int H, double s) : W_(W), H_(H), s_(s) {}

// ---------------------------------------------------------------------------
// computeFoldSpec — symmetric fold (§2)
//
// Symmetric fold over Z (crease parallel to X, FoldDirection::AlongX):
//   z_c = (H-1)s/2 ; a = (0, h, z_c) ; dHat = (1,0,0) ; nHat = (0,0,1)
//   mirror: (x, h + eps, 2 z_c - z)
// Symmetric fold over X (AlongZ): swap roles,
//   x_c = (W-1)s/2 ; nHat = (1,0,0) ; mirror x -> 2 x_c - x.
// Convention: the half with the LARGER coordinate moves (nHat points +).
// Crease-classified vertices (|r| <= tol) go into pinnedSet — S is pinned
// wholesale because pinning substitutes for table friction here (E9).
// ---------------------------------------------------------------------------
FoldSpec FoldingOptimizer::computeFoldSpec(double h, FoldDirection dir, double delta,
                                           double eps) const {
  (void)delta;   // arc-only parameter (§5); never influences the scoring targets (§2)
  FoldSpec spec;
  if (dir == FoldDirection::AlongX) {
    const double zc = (H_ - 1) * s_ * 0.5;
    spec.creasePoint  = Eigen::Vector3d(0.0, h, zc);
    spec.creaseDir    = Eigen::Vector3d(1.0, 0.0, 0.0);
    spec.creaseNormal = Eigen::Vector3d(0.0, 0.0, 1.0);
  } else {
    const double xc = (W_ - 1) * s_ * 0.5;
    spec.creasePoint  = Eigen::Vector3d(xc, h, 0.0);
    spec.creaseDir    = Eigen::Vector3d(0.0, 0.0, 1.0);
    spec.creaseNormal = Eigen::Vector3d(1.0, 0.0, 0.0);
  }
  classifyAndMirror(W_, H_, s_, h, eps, spec);

#ifndef NDEBUG
  // §4 paper checks, N = 10 symmetric over Z: the far five rows move
  // (indices 50..99); mirror of the far corner 99 = (0.9, h, 0.9) lands on
  // vertex 9 = (0.9, h, 0.0) lifted by eps.
  if (W_ == 10 && H_ == 10 && dir == FoldDirection::AlongX) {
    assert(spec.moving.size() == 50);
    for (size_t v : spec.moving) assert(v >= 50);
    auto it = std::find(spec.moving.begin(), spec.moving.end(), size_t(99));
    assert(it != spec.moving.end());
    const Eigen::Vector3d expect(0.9, h + eps, 0.0);
    assert((spec.mirrorTargets[size_t(it - spec.moving.begin())] - expect).norm() < 1e-9);
  }
#endif
  return spec;
}

// ---------------------------------------------------------------------------
// computeFoldSpec — diagonal fold (§2, general form)
//
// Given opposite corners startCorner (moves) and endCorner, the crease is the
// perpendicular bisector of the segment between them — for a square, the
// OTHER diagonal. a = cloth center = ((W-1)s/2, h, (H-1)s/2);
// nHat = normalize(pos(startCorner) - a) with the y component zeroed;
// dHat = normalize(yHat x nHat). Mirror: p' = p - 2 r_p nHat + eps*yHat.
// Canonical case (fold corner (N-1,N-1) onto (0,0)): crease is x + z = L,
// i.e. the lattice line i + j = N-1 (its vertices are pinned);
// M = { i + j > N-1 }; mirror(x, h, z) = (L - z, h + eps, L - x).
// ---------------------------------------------------------------------------
FoldSpec FoldingOptimizer::computeFoldSpec(double h, Corner startCorner, Corner endCorner,
                                           double delta, double eps) const {
  (void)delta;   // arc-only parameter (§5); never influences the scoring targets (§2)
  FoldSpec spec;
  if (W_ != H_) {
    std::cerr << "[FoldingOptimizer] diagonal fold requires a square grid (W == H); got "
              << W_ << "x" << H_ << " — returning empty spec." << std::endl;
    assert(W_ == H_ && "diagonal fold requires a square grid");
    return spec;
  }
  if (!cornersAreOpposite(startCorner, endCorner)) {
    std::cerr << "[FoldingOptimizer] diagonal fold requires OPPOSITE corners — returning empty spec."
              << std::endl;
    assert(cornersAreOpposite(startCorner, endCorner));
    return spec;
  }

  spec.creasePoint = Eigen::Vector3d((W_ - 1) * s_ * 0.5, h, (H_ - 1) * s_ * 0.5);
  Eigen::Vector3d n = cornerPosition(startCorner, W_, H_, s_, h) - spec.creasePoint;
  n.y() = 0.0;                       // in-plane normal toward the moving corner
  spec.creaseNormal = n.normalized();
  spec.creaseDir = kYHat.cross(spec.creaseNormal).normalized();
  classifyAndMirror(W_, H_, s_, h, eps, spec);

#ifndef NDEBUG
  // §4 paper checks, canonical N = 10 diagonal (FarRight -> NearLeft):
  // M = { i + j > 9 } has 9+8+...+1 = 45 vertices; the crease lattice line
  // i + j = 9 (indices 9, 18, ..., 90) is pinned; mirror of 99 = (0, h+eps, 0).
  if (W_ == 10 && H_ == 10 && startCorner == Corner::FarRight && endCorner == Corner::NearLeft) {
    assert(spec.moving.size() == 45);
    for (size_t k = 9; k <= 90; k += 9)
      assert(std::find(spec.pinnedSet.begin(), spec.pinnedSet.end(), k) != spec.pinnedSet.end());
    auto it = std::find(spec.moving.begin(), spec.moving.end(), size_t(99));
    assert(it != spec.moving.end());
    const Eigen::Vector3d expect(0.0, h + eps, 0.0);
    assert((spec.mirrorTargets[size_t(it - spec.moving.begin())] - expect).norm() < 1e-9);
  }
#endif
  return spec;
}

// ---------------------------------------------------------------------------
// gripSet — which moving vertices the "gripper" drives (§3, §4)
//
//   All        -> M itself.
//   FarEdge    -> symmetric: the outermost moving row/column (W grips);
//                 diagonal: the two outer edges of the moving triangle
//                 (canonical case: { v in M : i == H-1 or j == W-1 }).
//   TwoCorners -> symmetric: the two endpoints of that outer row/column;
//                 diagonal: not defined -> note + fall back to OneCorner.
//   OneCorner  -> the single far corner of M (max r, tie-break max u).
// Non-gripped moving vertices are left free (neither pinned nor controlled);
// that is the point of the ablation study. For the far-edge grip the
// semicircle is exactly the taut-string path (§9): the gripped edge stays at
// constant distance r_max equal to the material length, so the un-gripped
// interior is held on the straight taut column at its rigid position.
// ---------------------------------------------------------------------------
std::vector<size_t> FoldingOptimizer::gripSet(const FoldSpec& spec, GripMode mode) const {
  if (mode == GripMode::All || spec.moving.empty()) return spec.moving;

  const Eigen::Vector3d& n = spec.creaseNormal;
  const bool diagonal = std::abs(n.x()) > kTol && std::abs(n.z()) > kTol;
  // Extreme lattice row/column on the moving side (nHat points toward M).
  const int iOut = n.z() > kTol ? H_ - 1 : 0;
  const int jOut = n.x() > kTol ? W_ - 1 : 0;

  if (mode == GripMode::FarEdge) {
    std::vector<size_t> grips;
    for (size_t v : spec.moving) {
      const int i = static_cast<int>(v) / W_;
      const int j = static_cast<int>(v) % W_;
      const bool onOuter = diagonal ? (i == iOut || j == jOut)
                                    : (std::abs(n.z()) > kTol ? i == iOut : j == jOut);
      if (onOuter) grips.push_back(v);
    }
    return grips;
  }

  if (mode == GripMode::TwoCorners) {
    if (diagonal) {
      std::cerr << "[FoldingOptimizer] TwoCorners is not defined for the diagonal fold; "
                   "falling back to OneCorner." << std::endl;
      // fall through to OneCorner below
    } else {
      std::vector<size_t> grips;
      if (std::abs(n.z()) > kTol) {   // outer moving row -> its two endpoint corners
        grips.push_back(static_cast<size_t>(iOut) * static_cast<size_t>(W_));
        grips.push_back(static_cast<size_t>(iOut) * static_cast<size_t>(W_) + (W_ - 1));
      } else {                        // outer moving column -> its two endpoint corners
        grips.push_back(static_cast<size_t>(jOut));
        grips.push_back(static_cast<size_t>(H_ - 1) * static_cast<size_t>(W_) + jOut);
      }
      return grips;
    }
  }

  // OneCorner: the moving vertex farthest from the crease (max r); ties along
  // an outer row/column are broken by max u so the pick is deterministic.
  size_t best = spec.moving.front();
  double bestR = -1e300, bestU = -1e300;
  for (size_t v : spec.moving) {
    const int i = static_cast<int>(v) / W_;
    const int j = static_cast<int>(v) % W_;
    const Eigen::Vector3d p0(j * s_, spec.creasePoint.y(), i * s_);
    const double r = (p0 - spec.creasePoint).dot(spec.creaseNormal);
    const double u = (p0 - spec.creasePoint).dot(spec.creaseDir);
    if (r > bestR + kTol) { bestR = r; bestU = u; best = v; }
    else if (r > bestR - kTol && u > bestU) { bestU = u; best = v; }
  }
  return { best };
}

// ---------------------------------------------------------------------------
// arc — elliptical sweep about the crease with cosine easing (§5, §6;
//        FOLDING_FINAL_FORMULATION §elliptical-family)
//
// Each gripped vertex sweeps an arc about the crease. With
// c_v = p0 - r_v nHat (foot of the perpendicular on the crease; with crease
// slack delta: r_eff = r_v - delta and c_eff = p0 - r_eff nHat):
//
//   gamma(theta; alpha) = c_eff + r_eff cos(theta) nHat
//                         + alpha r_eff sin(theta) yHat
//                         + (theta/pi) eps yHat ,      theta: 0 -> pi
//
// alpha multiplies the VERTICAL (yHat) term ONLY — never the cos(theta) nHat
// term. Reason: the start (theta=0) and mirror-end (theta=pi) are fixed
// boundary conditions at horizontal offsets +r and -r; scaling cos(theta)
// would move those endpoints and the fold would not complete. So the
// horizontal semi-axis is pinned at r; only the height is free. alpha = 1 is
// the exact rigid semicircle (zero stretch); alpha < 1 lowers the path so the
// midpoint distance alpha*r < r leaves the material column slack (1-alpha)*r
// (it buckles: less lift, some stretch); alpha > 1 stretches the sheet
// (span > material length) and is Pareto-dominated. delta and alpha are
// independent: delta shrinks the radius (crease slack), alpha scales height.
//
// Timing (§6): cosine easing theta_k = (pi/2)(1 - cos(pi k/K)) at
// t_k = T k/K, k = 0..K — zero angular velocity at both endpoints.
//
// Why this path (§7, §9): a rigid rotation keeps every interior spring at
// rest length (zero stretch — unbeatable); crossing springs compress
// monotonically to a target-forced length, so this path minimizes peak
// elastic energy. Every moving point stays strictly above the pinned plane
// whenever it is horizontally over it (height h + r sin(theta) +
// (theta/pi) eps) — the ONLY interpenetration protection in an engine with
// no self-collision (E8). For the far-edge grip the semicircle is exactly
// the taut-string path (§9): the gripped edge stays at constant distance
// r_max equal to the material length, so the un-gripped interior is held on
// the straight taut column at its rigid position.
// ---------------------------------------------------------------------------
std::vector<Eigen::Vector3d> FoldingOptimizer::arc(const Eigen::Vector3d& p0, const FoldSpec& spec,
                                                   const FoldParams& params,
                                                   std::vector<double>& timesOut) const {
  timesOut.clear();
  std::vector<Eigen::Vector3d> pts;

  const double r = (p0 - spec.creasePoint).dot(spec.creaseNormal);
  const double rEff = r - params.delta;             // §5: crease slack shrinks the radius
  if (rEff <= kTol) {                               // degenerate: hold in place
    pts.push_back(p0);
    timesOut.push_back(0.0);
    return pts;
  }

  const Eigen::Vector3d cEff = p0 - rEff * spec.creaseNormal;
  for (int k = 0; k <= params.K; ++k) {
    const double theta = (kPi / 2.0) * (1.0 - std::cos(kPi * k / params.K));
    const double t = params.T * k / params.K;       // t_0 = 0 exactly (E6 trap 2)
    // alpha on yHat only (fixed endpoints pin the horizontal semi-axis at r).
    pts.push_back(cEff
                  + rEff * std::cos(theta) * spec.creaseNormal
                  + params.alpha * rEff * std::sin(theta) * kYHat
                  + (theta / kPi) * params.eps * kYHat);
    timesOut.push_back(t);
  }

  assert(timesOut.front() == 0.0);                  // E6: pre-start queries return the LAST
  for (size_t k = 1; k < timesOut.size(); ++k)      // waypoint, so times[0] must be exactly 0.0
    assert(timesOut[k] > timesOut[k - 1]);          // and times strictly increasing (NaN guard)
  return pts;
}

namespace {

// Shared tail of both optimize* entry points: grip -> arcs -> flat sorted list.
std::vector<TrajectoryPoint> buildTrajectory(const FoldingOptimizer& opt, const ClothMesh& cloth,
                                             const FoldSpec& spec, const FoldParams& params) {
  std::vector<TrajectoryPoint> out;
  for (size_t v : opt.gripSet(spec, params.grip)) {
    std::vector<double> times;
    const std::vector<Eigen::Vector3d> pts = opt.arc(cloth.getVertex(v).position, spec, params, times);
    for (size_t k = 0; k < pts.size(); ++k)
      out.push_back({ v, pts[k], times[k] });
  }
  std::sort(out.begin(), out.end(), [](const TrajectoryPoint& a, const TrajectoryPoint& b) {
    if (a.vertexIndex != b.vertexIndex) return a.vertexIndex < b.vertexIndex;
    return a.time < b.time;
  });
  return out;
}

} // namespace

// optimizeSymmetricFold / optimizeDiagonalFold (§3): compute the spec, take
// the grip set, generate arcs, and return one flat vector sorted by
// (vertexIndex, time). They do NOT touch the controller. The working height
// h is read off the cloth itself (drivers reposition the whole sheet to h
// after createGrid, which hardcodes y = 1.0).
std::vector<TrajectoryPoint> FoldingOptimizer::optimizeSymmetricFold(const ClothMesh& cloth,
                                                                     FoldDirection dir,
                                                                     double foldDuration,
                                                                     FoldParams params) {
  params.T = foldDuration;
  const double h = cloth.getVertex(0).position.y();
  const FoldSpec spec = computeFoldSpec(h, dir, params.delta, params.eps);
  return buildTrajectory(*this, cloth, spec, params);
}

std::vector<TrajectoryPoint> FoldingOptimizer::optimizeDiagonalFold(const ClothMesh& cloth,
                                                                    Corner startCorner,
                                                                    Corner endCorner,
                                                                    double foldDuration,
                                                                    FoldParams params) {
  params.T = foldDuration;
  const double h = cloth.getVertex(0).position.y();
  const FoldSpec spec = computeFoldSpec(h, startCorner, endCorner, params.delta, params.eps);
  return buildTrajectory(*this, cloth, spec, params);
}

// ---------------------------------------------------------------------------
// applyFoldingTrajectory — the E6 choke point (§3, §8)
//
// This must be the ONLY place in the codebase that calls setTrajectory,
// because setTrajectory carries two audited traps (E6):
//   Trap 1 (silent gain reset): setTrajectory internally calls
//     addPositionControl(v, positions[0]) with DEFAULT gain 1000 /
//     maxForce 100, silently overwriting anything configured earlier.
//     Fix: call addPositionControl(v, positions[0], gain, maxForce)
//     immediately AFTER setTrajectory — the trajectory only rewrites
//     targetPosition each frame, so gains installed afterwards stick.
//   Trap 2 (pre-start behavior): interpolateTrajectory returns the LAST
//     waypoint for any query time earlier than times[0]; the arcs therefore
//     guarantee times[0] == 0.0 exactly.
// Controller numbers (§8): stability requires gain < 2(1+d)/(d dt^2) ≈
// 42,000 at d = 0.9, dt = 0.01; tracking lag on a reference of speed V is
// e_ss = V(1-d)/(d gain dt) ≈ 2.5 mm at gain 2000 for the symmetric fold at
// T = 5. The defaults (gain 2000, maxForce 400) sit well inside that
// stability envelope.
// ---------------------------------------------------------------------------
void FoldingOptimizer::applyFoldingTrajectory(ClothController& controller,
                                              const std::vector<TrajectoryPoint>& trajectory,
                                              double gain, double maxForce) {
  // Group the flat vector by vertexIndex, preserving time order within each
  // vertex (optimize* outputs are (vertexIndex, time)-sorted; the map keeps
  // this correct even for unsorted callers).
  std::map<size_t, std::pair<std::vector<Eigen::Vector3d>, std::vector<double>>> perVertex;
  for (const TrajectoryPoint& tp : trajectory) {
    perVertex[tp.vertexIndex].first.push_back(tp.position);
    perVertex[tp.vertexIndex].second.push_back(tp.time);
  }
  for (const auto& kv : perVertex) {
    const size_t v = kv.first;
    const std::vector<Eigen::Vector3d>& positions = kv.second.first;
    const std::vector<double>& times = kv.second.second;
    controller.setTrajectory(v, positions, times, false);
    controller.addPositionControl(v, positions[0], gain, maxForce);  // re-arm gains (E6 trap 1)
  }
}

// ---------------------------------------------------------------------------
// gFoldTrajectory — gravity-based g-fold (van den Berg / Miller et al. eq. (1);
//                   FOLDING_FINAL_FORMULATION §g-fold)
//
// The real-world folding path: instead of a rigid arc, the gripper lifts the
// gripped corner along a 45-degree segment to peak height r above the crease,
// then descends mirrored — a triangle |perp| + height = r in the
// (fold-perpendicular, height) plane — while the UNGRIPPED cloth is left to
// hang under gravity (so this needs gravityScale > 1 and substeps to read).
//
// Eq. (1): a gripped point at fold-line distance y_v follows x = x_v (along
// the crease, unchanged), z = y_b (fold-perpendicular, swept y_v -> -y_v),
// y = y_v - |y_b| (height, peaking at y_v). Mapping to our frame with foot of
// perpendicular c_v = p0 - r nHat (on the crease at height h) and sweeping
// the perpendicular as y_b = r cos(theta):
//
//   gamma(theta) = c_v + r cos(theta) nHat
//                  + r (1 - |cos(theta)|) yHat        (triangular lift)
//                  + (theta/pi) eps yHat              (stacking gap)
//
// Endpoints match the arc's: theta=0 -> p0, theta=pi -> mirror + eps yHat, so
// the fold still LANDS at the same mirror target. This is a closed-form path;
// the papers' "optimization" is perception-side, not the trajectory. Times
// use the same cosine easing as arc (§6). Drives the OneCorner diagonal grip.
// ---------------------------------------------------------------------------
std::vector<TrajectoryPoint> FoldingOptimizer::gFoldTrajectory(const ClothMesh& cloth,
                                                               Corner startCorner,
                                                               Corner endCorner,
                                                               FoldParams params) {
  const double h = cloth.getVertex(0).position.y();
  const FoldSpec spec = computeFoldSpec(h, startCorner, endCorner, params.delta, params.eps);
  std::vector<TrajectoryPoint> out;
  if (spec.moving.empty()) return out;

  for (size_t v : gripSet(spec, params.grip)) {
    const Eigen::Vector3d p0 = cloth.getVertex(v).position;
    const double r = (p0 - spec.creasePoint).dot(spec.creaseNormal);
    const Eigen::Vector3d cFoot = p0 - r * spec.creaseNormal;   // on the crease, height h
    if (r <= kTol) { out.push_back({ v, p0, 0.0 }); continue; }  // degenerate: hold in place
    for (int k = 0; k <= params.K; ++k) {
      const double theta = (kPi / 2.0) * (1.0 - std::cos(kPi * k / params.K));  // §6 easing
      const double t = params.T * k / params.K;                  // t_0 = 0 exactly (E6)
      const double yb = r * std::cos(theta);                     // perpendicular sweep r -> -r
      const Eigen::Vector3d pos = cFoot + yb * spec.creaseNormal
                                  + r * (1.0 - std::abs(std::cos(theta))) * kYHat
                                  + (theta / kPi) * params.eps * kYHat;
      out.push_back({ v, pos, t });
    }
  }
  std::sort(out.begin(), out.end(), [](const TrajectoryPoint& a, const TrajectoryPoint& b) {
    if (a.vertexIndex != b.vertexIndex) return a.vertexIndex < b.vertexIndex;
    return a.time < b.time;
  });
  return out;
}

} // namespace ClothOpt

// folding.cpp — driver for the FoldingOptimizer (FOLDING_SOLUTION_v2 §10).
//
// Modes selected by argv:
//   (default)  GUI: Polyscope viewer modeled line-by-line on the working
//              idioms of optimization.cpp (init, YUp, registerSurfaceMesh,
//              state::userCallback, updateVertexPositions, show).
//   --eval     headless: run one complete fold with NO Polyscope call on the
//              code path, print one CSV line (plus header) to stdout, exit.
//   --grid     exhaustive (alpha, T) grid for the symmetric fold (CHANGE 7):
//              certifies the analytic arc optimum alpha* ~= 1.
//   --sweep    headless optimizer (CHANGE 8): 5-scalar CEM execution
//              calibration for symmetric folds, or spline + full-covariance
//              CEM over a 10-D gripper path for the diagonal one-corner fold.
//   --gfold    single headless run driving the OneCorner diagonal with the
//              gravity-based g-fold triangular path (CHANGE 10).
//
// The evaluation core runFold(...) is self-contained (its own ClothMesh,
// ClothController, integrator; no globals) so every mode reuses it.

#include "cloth.h"
#include "controller.h"
#include "integrator.h"
#include "folding_optimizer.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <random>
#include <string>
#include <vector>

#include <polyscope/polyscope.h>
#include <polyscope/surface_mesh.h>

using namespace ClothOpt;

// ===========================================================================
// Polyscope-free evaluation core
// ===========================================================================

enum class FoldType { SymZ, SymX, Diag };

struct RunConfig {
  FoldType fold = FoldType::SymZ;
  int N = 10;
  FoldParams P;          // grip lives inside FoldParams
  bool trace = false;
  bool gfold = false;    // CHANGE 10: drive the OneCorner diagonal with the g-fold path
  // CHANGE 8: an externally supplied trajectory (spline sweep). When non-null,
  // runFold drives the grip with this instead of generating an arc/g-fold.
  const std::vector<TrajectoryPoint>* customTraj = nullptr;
};

static constexpr double kSpacing = 0.1;   // §7 defaults: s = 0.1
static constexpr double kTol = 1e-9;

// gravityScale resolution (CHANGE 1): the boolean realGravity is kept for
// backward compatibility. If realGravity is on and gravityScale is still the
// default 1.0, reproduce the old E2b behavior (g x N*N); otherwise gravityScale
// is the continuous knob and is used directly.
static double resolveGravityScale(const FoldParams& P, int N) {
  if (P.realGravity && P.gravityScale == 1.0) return double(N) * double(N);
  return P.gravityScale;
}

// Linear interpolation of a per-vertex waypoint list, replicating the
// controller's interpolateTrajectory (E6): holds at positions.back() for
// t >= times.back(); with times[0] == 0 and queried t >= 0 the pre-start
// branch is never taken. Used to reconstruct the control target each step for
// the gripper-energy metric (CHANGE 4a).
static Eigen::Vector3d interpWaypoints(const std::vector<Eigen::Vector3d>& pos,
                                       const std::vector<double>& times, double t) {
  if (pos.empty()) return Eigen::Vector3d::Zero();
  if (t <= times.front()) return pos.front();
  if (t >= times.back()) return pos.back();
  size_t k = 1;
  while (k < times.size() && times[k] < t) ++k;
  const double span = times[k] - times[k - 1];
  const double a = span > 0 ? (t - times[k - 1]) / span : 0.0;
  return pos[k - 1] + a * (pos[k] - pos[k - 1]);
}

// One classified spring: crossing iff its endpoints have r of strictly
// opposite sign (both beyond tol) — classified ONCE at spec time (§10).
// Crossing springs are the ones forced away from rest by the fold itself
// (predicted crossing length law: s*cos(theta/2)); all others ("interior")
// should stay at rest length under the rigid-rotation trajectory (§7).
// The engine's distanceConstraints are structural only (rest s); there are no
// skip-one constraints — bendingConstraints (rest 2s) are the only other list.
// touching: exactly one endpoint has |r| <= kTol (crease-touching spring).
struct SpringRef { int v0, v1; double rest; bool crossing; bool touching; };

static std::vector<SpringRef> classifySprings(const ClothMesh& cloth, const FoldSpec& spec) {
  std::vector<SpringRef> out;
  auto rOf = [&](int v) {
    return (cloth.getVertex(static_cast<size_t>(v)).position - spec.creasePoint)
        .dot(spec.creaseNormal);
  };
  auto classify = [&](const std::vector<Edge>& edges) {
    for (const Edge& e : edges) {
      const double r0 = rOf(e.v0), r1 = rOf(e.v1);
      const bool crossing = (r0 > kTol && r1 < -kTol) || (r0 < -kTol && r1 > kTol);
      const bool touching = (std::abs(r0) <= kTol) != (std::abs(r1) <= kTol);
      out.push_back({ e.v0, e.v1, e.restLength, crossing, touching });
    }
  };
  classify(cloth.distanceConstraints);   // structural (rest s) — no skip-one constraints
  classify(cloth.bendingConstraints);    // bending (rest 2s) — E3: no shear springs exist
  return out;
}

static FoldSpec makeSpec(const RunConfig& cfg, const FoldingOptimizer& opt) {
  if (cfg.fold == FoldType::Diag)
    return opt.computeFoldSpec(cfg.P.h, Corner::FarRight, Corner::NearLeft, cfg.P.delta, cfg.P.eps);
  const FoldDirection dir = cfg.fold == FoldType::SymZ ? FoldDirection::AlongX
                                                       : FoldDirection::AlongZ;
  return opt.computeFoldSpec(cfg.P.h, dir, cfg.P.delta, cfg.P.eps);
}

// ---------------------------------------------------------------------------
// runFold — one complete headless rollout + metrics (§10)
//
// Rollout = T/dt fold steps + settle/dt settle steps, demo loop order
// (applyControls -> integrator step -> soft ground clamp -> [layer clamp]).
// Metrics:
//   M1 foldErrMean/Max : over ALL of M (not just grips),
//                        ||x_v(end) - mirrorTarget_v|| / L, L = (N-1)s.
//   M2 interiorStrainPeak : max over rollout & non-crossing springs of
//                        |len - rest| / rest — should be ~0 for a rigid
//                        rotation (§7); growth measures trajectory failure.
//   crossingMinLenRatio : min over rollout of crossing len / rest —
//                        companion of the predicted s*cos(theta/2) law (§9).
//   M6 interpenetrations : count of (moving vertex, frame) with
//                        (pos - a).nHat < -tol and pos.y < h - 0.002 (the
//                        flap dipping INTO the pinned layer, which nothing
//                        physical prevents, E8).
//   M8 residualKE : mean of 0.5|v|^2 over M during the final 1 s — end-state
//                        jitter (PBD vs controller fight, e.g. even-N crease).
//   E gripperEnergy : sum over steps and GRIPPED vertices of max(0, F.dx) —
//                        non-regenerative gripper work (raw; normalized to the
//                        arc baseline inside the cost lambdas, CHANGE 4a/7/8).
//   finite : false if any position goes non-finite (rollout aborts at once).
// Substeps (CHANGE 4b): the integrator is subdivided into P.substeps steps of
// dt/substeps per control step for stability at high stiffness / restored
// gravity. Gravity uses gravityScale (CHANGE 1/4c).
// ---------------------------------------------------------------------------
static FoldMetrics runFold(const RunConfig& cfg) {
  const int N = cfg.N;
  const FoldParams& P = cfg.P;
  FoldMetrics m;

  ClothMesh cloth;
  cloth.createGrid(N, N, kSpacing);
  // createGrid hardcodes y = 1.0 — reposition the sheet to working height h.
  for (size_t i = 0; i < cloth.getVertexCount(); ++i) {
    Eigen::Vector3d p = cloth.getVertex(i).position;
    p.y() = P.h;
    cloth.setVertexPosition(i, p);
  }
  cloth.properties.stiffness = P.stiffness;
  cloth.properties.bendingStiffness = P.bending;
  cloth.properties.damping = P.damping;
  // E1: the struct default gravity (0,0,-9.81) is HORIZONTAL vs the y=0
  // ground — must be set explicitly. E2/E2b: the integrator applies g/Nv;
  // gravityScale (resolved per CHANGE 1) restores true gravity (x N*N) or acts
  // as a continuous drape knob.
  const double gScale = resolveGravityScale(P, N);
  cloth.properties.gravity = gScale * Eigen::Vector3d(0.0, -9.81, 0.0);

  FoldingOptimizer opt(N, N, kSpacing);
  const FoldSpec spec = makeSpec(cfg, opt);
  if (spec.moving.empty()) { m.finite = false; return m; }

  const std::vector<SpringRef> springs = classifySprings(cloth, spec);  // at spec time, pre-motion

  std::vector<TrajectoryPoint> traj;
  {
    FoldParams params = P;
    if (cfg.gfold && cfg.fold != FoldType::Diag) {
      std::cerr << "[runFold] --gfold requires the diagonal fold" << std::endl;
      m.finite = false;
      return m;
    }
    if (cfg.customTraj)                       // CHANGE 8: spline sweep supplies its own path
      traj = *cfg.customTraj;
    else if (cfg.gfold)                       // CHANGE 10: g-fold triangular path (OneCorner diag)
      traj = opt.gFoldTrajectory(cloth, Corner::FarRight, Corner::NearLeft, params);
    else if (cfg.fold == FoldType::Diag)
      traj = opt.optimizeDiagonalFold(cloth, Corner::FarRight, Corner::NearLeft, P.T, params);
    else
      traj = opt.optimizeSymmetricFold(cloth, cfg.fold == FoldType::SymZ ? FoldDirection::AlongX
                                                                         : FoldDirection::AlongZ,
                                       P.T, params);
  }

  for (size_t v : spec.pinnedSet) cloth.pinVertex(static_cast<int>(v));

  // E10: the controlled set and the pinned set must be disjoint — a
  // controlled pinned vertex accumulates velocity forever, position frozen.
  const std::vector<size_t> grips = opt.gripSet(spec, P.grip);
  for (size_t g : grips) {
    if (cloth.pinned[g]) {
      std::cerr << "[runFold] FATAL: grip vertex " << g << " is pinned (E10)." << std::endl;
      std::abort();
    }
  }

  ClothController controller;
  auto integrator = std::make_unique<SemiImplicitEulerIntegrator>();
  opt.applyFoldingTrajectory(controller, traj, P.gain, P.maxForce);

  // CHANGE 4a: per-gripper waypoint lists, so we can reconstruct the control
  // target (and hence the clamped control force F) each step for the energy
  // metric. Only the gripped set contributes to gripperEnergy.
  std::map<size_t, std::pair<std::vector<Eigen::Vector3d>, std::vector<double>>> gripWaypoints;
  for (const TrajectoryPoint& tp : traj) {
    gripWaypoints[tp.vertexIndex].first.push_back(tp.position);
    gripWaypoints[tp.vertexIndex].second.push_back(tp.time);
  }

  const int foldSteps = static_cast<int>(std::lround(P.T / P.dt));
  const int settleSteps = static_cast<int>(std::lround(P.settle / P.dt));
  const int totalSteps = foldSteps + settleSteps;
  const int keWindow = std::min(settleSteps, static_cast<int>(std::lround(1.0 / P.dt)));

  std::ofstream traceFile;
  if (cfg.trace) {
    traceFile.open("crossing_trace.csv");
    traceFile << "frame,time,minCrossingRatio,meanCrossingRatio\n";
  }

  double keSum = 0.0;
  long keCount = 0;
  // Argmax tracking for interiorStrainPeak (FIX 5c).
  int argmaxV0 = -1, argmaxV1 = -1;
  double argmaxRest = 0.0;
  bool argmaxTouching = false;
  int argmaxStep = -1;

  for (int step = 0; step < totalSteps; ++step) {
    // CHANGE 4a: capture each gripper's position at the start of the control
    // step (this is the position the controller reads when it computes F).
    std::map<size_t, Eigen::Vector3d> posBefore;
    for (const auto& kv : gripWaypoints)
      posBefore[kv.first] = cloth.getVertex(kv.first).position;

    // Demo loop order (ground rules): controls, then physics, then clamps.
    // CHANGE 4b: the controller runs ONCE per control step; the integrator is
    // subdivided into P.substeps sub-steps of dt/substeps (the Simulation
    // demo's own pattern, required for stability at high stiffness / restored
    // gravity).
    controller.applyControls(cloth, P.dt);
    const int substeps = std::max(1, P.substeps);
    for (int sub = 0; sub < substeps; ++sub) integrator->step(cloth, P.dt / substeps);

    // E8: soft ground clamp replicated from the Optimization demo's frame
    // loop — y < 0.01 snaps to 0.01 with vy *= 0.3 (the integrator's own
    // clamp at y < 0 has no velocity handling at all).
    for (size_t i = 0; i < cloth.getVertexCount(); ++i) {
      Eigen::Vector3d p = cloth.getVertex(i).position;
      if (p.y() < 0.01) {
        p.y() = 0.01;
        cloth.setVertexPosition(i, p);
        cloth.velocities[i].y() *= 0.3;
      }
    }

    // Virtual layer clamp (§10): stand-in for the missing cloth-cloth layer
    // contact, generalizing the engine's own ground-clamp idiom to the plane
    // of the pinned half — a moving vertex that has crossed the crease
    // ((pos - a).nHat < -tol) may not sink below the stacking height h + eps.
    if (P.layerClamp) {
      for (size_t v : spec.moving) {
        Eigen::Vector3d p = cloth.getVertex(v).position;
        const double rNow = (p - spec.creasePoint).dot(spec.creaseNormal);
        if (rNow < -kTol && p.y() < P.h + P.eps) {
          p.y() = P.h + P.eps;
          cloth.setVertexPosition(v, p);
          cloth.velocities[v].y() *= 0.3;
        }
      }
    }

    // CHANGE 4a: gripper energy E — non-regenerative work the controller did
    // this step, summed over gripped vertices. The controller applied
    // F = clamp(gain*(target - pos_before), maxForce) as v += F*dt (E5);
    // reconstruct F from the same target it used (its clock reads (step+1)*dt
    // after applyControls) and dx = pos_after - pos_before. Only positive
    // work counts (max(0, F.dx)): energy the gripper injects, not what the
    // cloth returns.
    const double controllerTime = (step + 1) * P.dt;
    for (const auto& kv : gripWaypoints) {
      const size_t v = kv.first;
      const Eigen::Vector3d target = interpWaypoints(kv.second.first, kv.second.second, controllerTime);
      Eigen::Vector3d F = P.gain * (target - posBefore[v]);
      if (F.norm() > P.maxForce) F = F.normalized() * P.maxForce;
      const Eigen::Vector3d dx = cloth.getVertex(v).position - posBefore[v];
      m.gripperEnergy += std::max(0.0, F.dot(dx));
    }

    // ---- metrics block (§10; definitions in the function header) ----
    for (size_t i = 0; i < cloth.getVertexCount(); ++i) {
      if (!cloth.getVertex(i).position.allFinite()) {
        m.finite = false;
        return m;                                     // abort the rollout immediately
      }
    }

    double minCross = std::numeric_limits<double>::infinity();
    double sumCross = 0.0;
    int nCross = 0;
    for (const SpringRef& sp : springs) {
      const double len = (cloth.getVertex(static_cast<size_t>(sp.v0)).position -
                          cloth.getVertex(static_cast<size_t>(sp.v1)).position).norm();
      if (sp.crossing) {
        const double ratio = len / sp.rest;
        minCross = std::min(minCross, ratio);
        sumCross += ratio;
        ++nCross;
        m.crossingMinLenRatio = std::min(m.crossingMinLenRatio, ratio);
      } else {
        const double strain = std::abs(len - sp.rest) / sp.rest;
        if (strain > m.interiorStrainPeak) {
          m.interiorStrainPeak = strain;
          argmaxV0 = sp.v0; argmaxV1 = sp.v1;
          argmaxRest = sp.rest; argmaxTouching = sp.touching;
          argmaxStep = step;
        }
        if (sp.touching) m.touchStrainPeak = std::max(m.touchStrainPeak, strain);
      }
    }

    for (size_t v : spec.moving) {
      const Eigen::Vector3d& p = cloth.getVertex(v).position;
      const double rNow = (p - spec.creasePoint).dot(spec.creaseNormal);
      if (rNow < -kTol && p.y() < P.h - 0.002) ++m.interpenetrations;
    }

    if (step >= totalSteps - keWindow) {
      for (size_t v : spec.moving) keSum += 0.5 * cloth.velocities[v].squaredNorm();
      keCount += static_cast<long>(spec.moving.size());
    }

    if (traceFile.is_open() && nCross > 0)
      traceFile << step << ',' << (step + 1) * P.dt << ',' << minCross << ','
                << sumCross / nCross << '\n';
  }

  m.residualKE = keCount > 0 ? keSum / static_cast<double>(keCount) : 0.0;

  // Print the spring that dominates interiorStrainPeak (FIX 5c).
  std::cerr << "[strain-argmax] v0=" << argmaxV0 << " v1=" << argmaxV1
            << " rest=" << argmaxRest << " touching=" << (argmaxTouching ? 1 : 0)
            << " step=" << argmaxStep << " peak=" << m.interiorStrainPeak << std::endl;

  // M1: end-state placement error over ALL of M, as a fraction of side L.
  const double L = (N - 1) * kSpacing;
  double sum = 0.0, mx = 0.0;
  for (size_t k = 0; k < spec.moving.size(); ++k) {
    const double err =
        (cloth.getVertex(spec.moving[k]).position - spec.mirrorTargets[k]).norm() / L;
    sum += err;
    mx = std::max(mx, err);
  }
  m.foldErrMean = sum / static_cast<double>(spec.moving.size());
  m.foldErrMax = mx;
  return m;
}

// ===========================================================================
// Startup self-test — the §4 paper checks for N = 10, run unconditionally
// (Release builds compile assert() out via -DNDEBUG). Output goes to stderr
// so --eval's CSV-on-stdout contract stays clean.
// ===========================================================================
static bool selfTest() {
  bool ok = true;
  const double h = 0.25, eps = 0.02;
  FoldingOptimizer opt(10, 10, kSpacing);

  const FoldSpec sym = opt.computeFoldSpec(h, FoldDirection::AlongX, 0.0, eps);
  ok &= sym.moving.size() == 50;
  for (size_t v : sym.moving) ok &= v >= 50;
  {
    auto it = std::find(sym.moving.begin(), sym.moving.end(), size_t(99));
    ok &= it != sym.moving.end();
    if (it != sym.moving.end()) {
      const Eigen::Vector3d expect(0.9, h + eps, 0.0);   // vertex 9 + (0, eps, 0)
      ok &= (sym.mirrorTargets[size_t(it - sym.moving.begin())] - expect).norm() < 1e-9;
    }
  }

  const FoldSpec diag = opt.computeFoldSpec(h, Corner::FarRight, Corner::NearLeft, 0.0, eps);
  ok &= diag.moving.size() == 45;
  for (size_t k = 9; k <= 90; k += 9)
    ok &= std::find(diag.pinnedSet.begin(), diag.pinnedSet.end(), k) != diag.pinnedSet.end();
  {
    auto it = std::find(diag.moving.begin(), diag.moving.end(), size_t(99));
    ok &= it != diag.moving.end();
    if (it != diag.moving.end()) {
      const Eigen::Vector3d expect(0.0, h + eps, 0.0);
      ok &= (diag.mirrorTargets[size_t(it - diag.moving.begin())] - expect).norm() < 1e-9;
    }
  }

  std::cerr << "[selftest] fold-spec invariants (N=10, sym+diag): "
            << (ok ? "OK" : "FAILED") << std::endl;
  return ok;
}

// ===========================================================================
// Flag parsing (shared by --eval and --sweep)
// ===========================================================================
struct ParsedArgs {
  RunConfig cfg;
  std::string foldName = "sym";
  std::string gripName = "all";
  bool tSet = false;
  bool gripSet = false;   // FIX 3: track whether --grip was explicitly passed
  bool foldSet = false;  // FIX 3: track whether --fold was explicitly passed
  // CHANGE 7: exhaustive symmetric grid controls (alpha x T) + energy weight.
  int nAlpha = 25, nT = 20;
  double wE = 0.2;
};

static bool parseArgs(int argc, char** argv, ParsedArgs& out) {
  auto need = [&](int i) -> const char* {
    if (i + 1 >= argc) {
      std::cerr << "missing value for " << argv[i] << std::endl;
      std::exit(2);
    }
    return argv[i + 1];
  };
  for (int i = 1; i < argc; ++i) {
    const std::string a = argv[i];
    if (a == "--eval" || a == "--sweep" || a == "--grid") continue;  // mode flags (no value)
    else if (a == "--gfold") { out.cfg.gfold = true; }               // CHANGE 10 (no value)
    else if (a == "--fold") { out.foldName = need(i); out.foldSet = true; ++i; }
    else if (a == "--grip") { out.gripName = need(i); out.gripSet = true; ++i; }
    else if (a == "--N") { out.cfg.N = std::atoi(need(i)); ++i; }
    else if (a == "--T") { out.cfg.P.T = std::atof(need(i)); out.tSet = true; ++i; }
    else if (a == "--gain") { out.cfg.P.gain = std::atof(need(i)); ++i; }
    else if (a == "--maxForce") { out.cfg.P.maxForce = std::atof(need(i)); ++i; }
    else if (a == "--eps") { out.cfg.P.eps = std::atof(need(i)); ++i; }
    else if (a == "--delta") { out.cfg.P.delta = std::atof(need(i)); ++i; }
    else if (a == "--stiff") { out.cfg.P.stiffness = std::atof(need(i)); ++i; }
    else if (a == "--bend") { out.cfg.P.bending = std::atof(need(i)); ++i; }
    else if (a == "--damp") { out.cfg.P.damping = std::atof(need(i)); ++i; }
    else if (a == "--realGravity") { out.cfg.P.realGravity = std::atoi(need(i)) != 0; ++i; }
    else if (a == "--layerClamp") { out.cfg.P.layerClamp = std::atoi(need(i)) != 0; ++i; }
    else if (a == "--settle") { out.cfg.P.settle = std::atof(need(i)); ++i; }
    else if (a == "--trace") { out.cfg.trace = std::atoi(need(i)) != 0; ++i; }
    // CHANGE 6: new continuous knobs.
    else if (a == "--alpha") { out.cfg.P.alpha = std::atof(need(i)); ++i; }
    else if (a == "--substeps") { out.cfg.P.substeps = std::atoi(need(i)); ++i; }
    else if (a == "--gravityScale") { out.cfg.P.gravityScale = std::atof(need(i)); ++i; }
    // CHANGE 7: grid controls.
    else if (a == "--Nalpha") { out.nAlpha = std::atoi(need(i)); ++i; }
    else if (a == "--NT") { out.nT = std::atoi(need(i)); ++i; }
    else if (a == "--wE") { out.wE = std::atof(need(i)); ++i; }
    else { std::cerr << "unknown flag: " << a << std::endl; return false; }
  }

  // FIX 3: --gfold forces the diagonal fold and (if --grip was not explicitly
  // given) the corner grip so the run is self-consistent.
  if (out.cfg.gfold) {
    if (!out.foldSet) out.foldName = "diag";
    if (!out.gripSet) out.gripName = "corner";
  }

  if (out.foldName == "sym") out.cfg.fold = FoldType::SymZ;
  else if (out.foldName == "symx") out.cfg.fold = FoldType::SymX;
  else if (out.foldName == "diag") out.cfg.fold = FoldType::Diag;
  else { std::cerr << "--fold must be sym|symx|diag" << std::endl; return false; }

  if (out.gripName == "all") out.cfg.P.grip = GripMode::All;
  else if (out.gripName == "edge") out.cfg.P.grip = GripMode::FarEdge;
  else if (out.gripName == "corners") out.cfg.P.grip = GripMode::TwoCorners;
  else if (out.gripName == "corner") out.cfg.P.grip = GripMode::OneCorner;
  else { std::cerr << "--grip must be all|edge|corners|corner" << std::endl; return false; }

  // §7: T defaults to 5 for symmetric, 6 for the diagonal.
  if (!out.tSet && out.cfg.fold == FoldType::Diag) out.cfg.P.T = 6.0;
  return true;
}

static void printCsv(const ParsedArgs& pa, const FoldMetrics& m) {
  const FoldParams& P = pa.cfg.P;
  // CHANGE 5: added alpha, substeps, gravityScale (the resolved effective value
  // that actually ran, per CHANGE 1) and gripperEnergy. Existing columns kept.
  std::cout << "fold,grip,N,T,gain,maxForce,eps,delta,stiff,bend,damp,realGravity,layerClamp,"
               "settle,alpha,substeps,gravityScale,foldErrMean,foldErrMax,interiorStrainPeak,"
               "crossingMinLenRatio,interpenetrations,residualKE,finite,gripperEnergy,touchStrainPeak\n";
  std::cout << std::setprecision(6)
            << pa.foldName << ',' << pa.gripName << ',' << pa.cfg.N << ',' << P.T << ','
            << P.gain << ',' << P.maxForce << ',' << P.eps << ',' << P.delta << ','
            << P.stiffness << ',' << P.bending << ',' << P.damping << ','
            << (P.realGravity ? 1 : 0) << ',' << (P.layerClamp ? 1 : 0) << ',' << P.settle << ','
            << P.alpha << ',' << P.substeps << ',' << resolveGravityScale(P, pa.cfg.N) << ','
            << m.foldErrMean << ',' << m.foldErrMax << ',' << m.interiorStrainPeak << ','
            << m.crossingMinLenRatio << ',' << m.interpenetrations << ',' << m.residualKE << ','
            << (m.finite ? 1 : 0) << ',' << m.gripperEnergy << ',' << m.touchStrainPeak << '\n';
}

static int runEval(const ParsedArgs& pa) {
  printCsv(pa, runFold(pa.cfg));
  return 0;
}

// ===========================================================================
// --grid : EXHAUSTIVE 2D grid over (alpha, T) for the SYMMETRIC fold
// (FOLDING_FINAL_FORMULATION §elliptical-family). Rationale for an exhaustive
// grid: the symmetric search is only 2-D, so a grid is exact and fully
// interpretable — it maps the whole Pareto surface and CERTIFIES that the
// optimum lands at alpha* ~= 1 (the analytic rigid arc) rather than assuming
// it. The objective is a black-box, non-smooth rollout functional (force
// clamp + PBD projection + penalty), so no usable gradient exists; rollouts
// are cheap (~65 ms), so dense sampling is practical.
//
// Cost:  J = 1.0*D + 1.0*S + wE*(E / E_ref) + 1000*(C>0 || !finite),
// where E_ref is the gripper energy of the (alpha=1, SAME T) run, so E is
// normalized to the rigid arc at that T. Prints one CSV row per cell and the
// best (min-J) cell at the end.
// ===========================================================================
static int runGrid(const ParsedArgs& base) {
  const int nA = std::max(2, base.nAlpha);
  const int nT = std::max(2, base.nT);
  const double wE = base.wE;

  std::cout << "alpha,T,D,S,E,C,resKE,finite,J\n";
  double bestJ = std::numeric_limits<double>::infinity();
  double bestAlpha = 1.0, bestT = base.cfg.P.T;

  for (int it = 0; it < nT; ++it) {
    const double T = 2.0 + (10.0 - 2.0) * it / (nT - 1);

    // E_ref: the rigid-arc (alpha=1) gripper energy at THIS T.
    RunConfig ref = base.cfg;
    ref.trace = false;
    ref.P.alpha = 1.0;
    ref.P.T = T;
    const FoldMetrics refM = runFold(ref);
    const double eRef = refM.gripperEnergy > 1e-12 ? refM.gripperEnergy : 1.0;

    for (int ia = 0; ia < nA; ++ia) {
      const double alpha = 0.5 + (1.2 - 0.5) * ia / (nA - 1);
      RunConfig c = base.cfg;
      c.trace = false;
      c.P.alpha = alpha;
      c.P.T = T;
      const FoldMetrics m = runFold(c);
      const double D = m.foldErrMean, S = m.interiorStrainPeak, E = m.gripperEnergy;
      const int C = m.interpenetrations;
      double J = 1.0 * D + 1.0 * S + wE * (E / eRef);
      if (C > 0 || !m.finite) J += 1000.0;
      std::cout << std::setprecision(6) << alpha << ',' << T << ',' << D << ',' << S << ','
                << E << ',' << C << ',' << m.residualKE << ',' << (m.finite ? 1 : 0) << ',' << J << '\n';
      if (J < bestJ) { bestJ = J; bestAlpha = alpha; bestT = T; }
    }
  }

  std::cout << "[grid] BEST J=" << std::setprecision(6) << bestJ
            << " at alpha=" << bestAlpha << " T=" << bestT << std::endl;
  return 0;
}

// ===========================================================================
// --grid --fold diag : exhaustive (alpha, T) grid for the DIAGONAL fold,
// one-corner (far-corner) grip. Mirrors runGrid but with FoldType::Diag and
// GripMode::OneCorner forced. The same alpha (elliptical arc height) and T
// (fold duration) parametrize the one-corner trajectory via the shared arc()
// function, so the 2-D grid has the same semantic as for the symmetric fold.
// Cost:  J = 1.0*D + 1.0*S + wE*(E / E_ref) + 1000*(C>0 || !finite)
// where E_ref = gripper energy of the (alpha=1, SAME T) diagonal run.
// ===========================================================================
static int runGridDiagonal(const ParsedArgs& base) {
  const int nA = std::max(2, base.nAlpha);
  const int nT = std::max(2, base.nT);
  const double wE = base.wE;

  std::cout << "alpha,T,D,S,E,C,resKE,finite,J\n";
  double bestJ = std::numeric_limits<double>::infinity();
  double bestAlpha = 1.0, bestT = base.cfg.P.T;

  for (int it = 0; it < nT; ++it) {
    const double T = 2.0 + (10.0 - 2.0) * it / (nT - 1);

    // E_ref: the rigid-arc (alpha=1, OneCorner diagonal) gripper energy at THIS T.
    RunConfig ref = base.cfg;
    ref.fold = FoldType::Diag;
    ref.P.grip = GripMode::OneCorner;
    ref.trace = false;
    ref.gfold = false;
    ref.customTraj = nullptr;
    ref.P.alpha = 1.0;
    ref.P.T = T;
    const FoldMetrics refM = runFold(ref);
    const double eRef = refM.gripperEnergy > 1e-12 ? refM.gripperEnergy : 1.0;

    for (int ia = 0; ia < nA; ++ia) {
      const double alpha = 0.5 + (1.2 - 0.5) * ia / (nA - 1);
      RunConfig c = base.cfg;
      c.fold = FoldType::Diag;
      c.P.grip = GripMode::OneCorner;
      c.trace = false;
      c.gfold = false;
      c.customTraj = nullptr;
      c.P.alpha = alpha;
      c.P.T = T;
      const FoldMetrics m = runFold(c);
      const double D = m.foldErrMean, S = m.interiorStrainPeak, E = m.gripperEnergy;
      const int C = m.interpenetrations;
      double J = 1.0 * D + 1.0 * S + wE * (E / eRef);
      if (C > 0 || !m.finite) J += 1000.0;
      std::cout << std::setprecision(6) << alpha << ',' << T << ',' << D << ',' << S << ','
                << E << ',' << C << ',' << m.residualKE << ',' << (m.finite ? 1 : 0) << ',' << J << '\n';
      if (J < bestJ) { bestJ = J; bestAlpha = alpha; bestT = T; }
    }
  }

  std::cout << "[grid] BEST J=" << std::setprecision(6) << bestJ
            << " at alpha=" << bestAlpha << " T=" << bestT << std::endl;
  return 0;
}

// ===========================================================================
// CHANGE 8 — diagonal one-corner optimizer: cubic Catmull-Rom gripper path +
// full-covariance CEM.
//
// Why genuine optimization here (not a grid): the one-corner diagonal has NO
// closed-form ideal path — a single gripper cannot hold a 2-D sheet rigid, so
// the trajectory that best trades fold error against stretch and energy must
// be searched. The space is ~10-D (9 control-point coords + T) with
// correlated parameters, where an exhaustive grid is hopeless; a
// derivative-free method with a learned full covariance is the right tool
// (the objective is a non-smooth black-box rollout, no gradient). A
// differentiable simulator (e.g. SoftMAC) would enable gradient-based
// trajectory optimization — future work, not this task.
// ===========================================================================

// Uniform Catmull-Rom segment P1->P2 with neighbors P0, P3 (endpoints
// duplicated at the ends). Standard basis form.
static Eigen::Vector3d catmullRom(const Eigen::Vector3d& p0, const Eigen::Vector3d& p1,
                                  const Eigen::Vector3d& p2, const Eigen::Vector3d& p3, double t) {
  const double t2 = t * t, t3 = t2 * t;
  return 0.5 * ((2.0 * p1) + (-p0 + p2) * t + (2.0 * p0 - 5.0 * p1 + 4.0 * p2 - p3) * t2 +
                (-p0 + 3.0 * p1 - 3.0 * p2 + p3) * t3);
}

// buildSplineTrajectory (CHANGE 8a): a cubic Catmull-Rom spline through
// {start, p1, p2, p3, end} (start = gripped corner's rest position, end = its
// mirror target, p1..p3 = 3 free interior control points). Discretized into
// K+1 waypoints whose PATH parameter is cosine-eased (u_k = (1-cos(pi k/K))/2)
// so speed vanishes at both ends, at strictly increasing times t_k = T k/K
// (t_0 = 0 exactly, E6). vertexIndex added to the spec'd signature so the flat
// TrajectoryPoint list carries its vertex.
static std::vector<TrajectoryPoint> buildSplineTrajectory(size_t vertexIndex,
                                                          const Eigen::Vector3d& cornerRest,
                                                          const Eigen::Vector3d& cornerMirror,
                                                          const Eigen::Vector3d& p1,
                                                          const Eigen::Vector3d& p2,
                                                          const Eigen::Vector3d& p3,
                                                          double T, int K) {
  const std::array<Eigen::Vector3d, 5> P = { cornerRest, p1, p2, p3, cornerMirror };
  std::vector<TrajectoryPoint> out;
  out.reserve(size_t(K) + 1);
  for (int k = 0; k <= K; ++k) {
    const double u = 0.5 * (1.0 - std::cos(M_PI * k / K));   // cosine-eased path parameter
    const double s = std::min(std::max(u, 0.0), 1.0) * 4.0;  // 4 segments
    int seg = std::min(3, static_cast<int>(std::floor(s)));
    const double lt = s - seg;                               // local t in [0,1]
    const Eigen::Vector3d& a = P[size_t(std::max(0, seg - 1))];
    const Eigen::Vector3d& b = P[size_t(seg)];
    const Eigen::Vector3d& c = P[size_t(seg + 1)];
    const Eigen::Vector3d& d = P[size_t(std::min(4, seg + 2))];
    out.push_back({ vertexIndex, catmullRom(a, b, c, d, lt), T * k / K });
  }
  return out;
}

static int runSweepDiagonal(const ParsedArgs& base) {
  const int N = base.cfg.N;
  const double h = base.cfg.P.h, s = kSpacing;
  const int K = base.cfg.P.K;

  FoldingOptimizer opt(N, N, s);
  const FoldSpec spec = opt.computeFoldSpec(h, Corner::FarRight, Corner::NearLeft,
                                            base.cfg.P.delta, base.cfg.P.eps);
  const std::vector<size_t> grips = opt.gripSet(spec, GripMode::OneCorner);
  if (grips.empty()) { std::cerr << "[sweep-diag] empty grip set" << std::endl; return 1; }
  const size_t gripV = grips[0];
  const int gi = static_cast<int>(gripV) / N, gj = static_cast<int>(gripV) % N;
  const Eigen::Vector3d cornerRest(gj * s, h, gi * s);
  Eigen::Vector3d cornerMirror = cornerRest;
  for (size_t k = 0; k < spec.moving.size(); ++k)
    if (spec.moving[k] == gripV) cornerMirror = spec.mirrorTargets[k];

  // Seed the 3 interior control points ON the analytic alpha=1 arc.
  FoldParams arcParams = base.cfg.P;
  arcParams.alpha = 1.0;
  arcParams.T = 6.0;
  std::vector<double> arcTimes;
  const std::vector<Eigen::Vector3d> arcPts = opt.arc(cornerRest, spec, arcParams, arcTimes);
  const std::array<Eigen::Vector3d, 3> seedP = { arcPts[size_t(K / 4)], arcPts[size_t(K / 2)],
                                                 arcPts[size_t(3 * K / 4)] };

  // E_ref: gripper energy of the analytic arc (OneCorner diagonal, T=6).
  auto arcRefConfig = [&]() {
    RunConfig c = base.cfg;
    c.fold = FoldType::Diag;
    c.P.grip = GripMode::OneCorner;
    c.P.alpha = 1.0;
    c.P.T = 6.0;
    c.gfold = false;
    c.customTraj = nullptr;
    return c;
  };
  const double eRef = std::max(1e-12, runFold(arcRefConfig()).gripperEnergy);

  // 10-D box: each control coord in seed +/- 0.45 m; T in [2,10]. The arc seed
  // maps to z = 0.5 in every dimension.
  constexpr int D = 10;
  std::array<double, D> lo{}, hi{};
  for (int p = 0; p < 3; ++p)
    for (int c = 0; c < 3; ++c) {
      lo[size_t(3 * p + c)] = seedP[size_t(p)][c] - 0.45;
      hi[size_t(3 * p + c)] = seedP[size_t(p)][c] + 0.45;
    }
  lo[9] = 2.0; hi[9] = 10.0;

  auto zToTraj = [&](const std::array<double, D>& z) {
    Eigen::Vector3d p1, p2, p3;
    for (int c = 0; c < 3; ++c) {
      p1[c] = lo[size_t(c)]     + z[size_t(c)]     * (hi[size_t(c)]     - lo[size_t(c)]);
      p2[c] = lo[size_t(3 + c)] + z[size_t(3 + c)] * (hi[size_t(3 + c)] - lo[size_t(3 + c)]);
      p3[c] = lo[size_t(6 + c)] + z[size_t(6 + c)] * (hi[size_t(6 + c)] - lo[size_t(6 + c)]);
    }
    const double T = lo[9] + z[9] * (hi[9] - lo[9]);
    return std::make_pair(buildSplineTrajectory(gripV, cornerRest, cornerMirror, p1, p2, p3, T, K), T);
  };
  auto cost = [&](const std::array<double, D>& z, FoldMetrics& m) {
    auto tr = zToTraj(z);
    RunConfig c = base.cfg;
    c.fold = FoldType::Diag;
    c.P.grip = GripMode::OneCorner;
    c.P.T = tr.second;
    c.gfold = false;
    c.customTraj = &tr.first;
    m = runFold(c);
    double J = m.foldErrMean + m.interiorStrainPeak + base.wE * (m.gripperEnergy / eRef);
    if (m.interpenetrations > 0 || !m.finite) J += 1000.0;
    return J;
  };

  std::mt19937 rng(1234);
  std::uniform_real_distribution<double> uni(0.0, 1.0);
  std::normal_distribution<double> gauss(0.0, 1.0);

  // Pre-scan: LHS(128) (seed included at z = 0.5) -> keep best 16 as the
  // initial elite pool.
  const int NS = 128;
  std::array<std::vector<int>, D> perms;
  for (auto& pm : perms) {
    pm.resize(size_t(NS));
    for (int k = 0; k < NS; ++k) pm[size_t(k)] = k;
    std::shuffle(pm.begin(), pm.end(), rng);
  }
  std::vector<std::pair<double, std::array<double, D>>> pool;
  {
    std::array<double, D> seedZ; seedZ.fill(0.5);   // the analytic arc
    FoldMetrics m; pool.emplace_back(cost(seedZ, m), seedZ);
  }
  for (int k = 0; k < NS; ++k) {
    std::array<double, D> z;
    for (int d = 0; d < D; ++d) z[size_t(d)] = (perms[size_t(d)][size_t(k)] + uni(rng)) / NS;
    FoldMetrics m; pool.emplace_back(cost(z, m), z);
  }
  std::sort(pool.begin(), pool.end(), [](const auto& a, const auto& b) { return a.first < b.first; });
  pool.resize(16);
  auto best = pool.front();
  std::cout << "[sweep-diag] LHS best J=" << std::setprecision(5) << best.first << std::endl;

  // Full-covariance CEM: fit a 10-D Gaussian to the elites, sample via its
  // Cholesky factor (with a small diagonal jitter for positive-definiteness),
  // refit on the elite 8 each iteration. (The spec permits full-covariance CEM
  // in place of CMA-ES; this IS that, labelled.)
  Eigen::VectorXd mean(D);
  Eigen::MatrixXd cov(D, D);
  auto fit = [&](const std::vector<std::pair<double, std::array<double, D>>>& elites) {
    mean.setZero();
    for (const auto& e : elites)
      for (int d = 0; d < D; ++d) mean[d] += e.second[size_t(d)];
    mean /= double(elites.size());
    cov.setZero();
    for (const auto& e : elites) {
      Eigen::VectorXd x(D);
      for (int d = 0; d < D; ++d) x[d] = e.second[size_t(d)];
      cov += (x - mean) * (x - mean).transpose();
    }
    cov /= double(elites.size());
    cov += Eigen::MatrixXd::Identity(D, D) * 1e-4;   // jitter: PD + sigma floor
  };
  fit(pool);

  for (int iter = 1; iter <= 12; ++iter) {
    Eigen::LLT<Eigen::MatrixXd> llt(cov);
    const Eigen::MatrixXd L = llt.matrixL();
    std::vector<std::pair<double, std::array<double, D>>> gen;
    for (int k = 0; k < 32; ++k) {
      Eigen::VectorXd n(D);
      for (int d = 0; d < D; ++d) n[d] = gauss(rng);
      const Eigen::VectorXd x = mean + L * n;
      std::array<double, D> z;
      for (int d = 0; d < D; ++d) z[size_t(d)] = std::clamp(x[d], 0.0, 1.0);
      FoldMetrics m; gen.emplace_back(cost(z, m), z);
    }
    std::sort(gen.begin(), gen.end(), [](const auto& a, const auto& b) { return a.first < b.first; });
    if (gen.front().first < best.first) best = gen.front();
    gen.resize(8);                                   // elite 8
    fit(gen);
    std::cout << "[sweep-diag] iter " << iter << " best J=" << std::setprecision(5)
              << best.first << std::endl;
  }

  // Report the best 10-parameter trajectory and its full metric line (FIX 1).
  auto tr = zToTraj(best.second);
  RunConfig c = base.cfg;
  c.fold = FoldType::Diag;
  c.P.grip = GripMode::OneCorner;
  c.P.T = tr.second;
  c.customTraj = &tr.first;

  // Decode winning control points from best.second using the same lo/hi mapping.
  Eigen::Vector3d winP1, winP2, winP3;
  for (int ci = 0; ci < 3; ++ci) {
    winP1[ci] = lo[size_t(ci)]     + best.second[size_t(ci)]     * (hi[size_t(ci)]     - lo[size_t(ci)]);
    winP2[ci] = lo[size_t(3 + ci)] + best.second[size_t(3 + ci)] * (hi[size_t(3 + ci)] - lo[size_t(3 + ci)]);
    winP3[ci] = lo[size_t(6 + ci)] + best.second[size_t(6 + ci)] * (hi[size_t(6 + ci)] - lo[size_t(6 + ci)]);
  }
  std::cout << std::setprecision(9)
            << "[sweep-diag] start   = (" << cornerRest.x()   << ", " << cornerRest.y()   << ", " << cornerRest.z()   << ")\n"
            << "[sweep-diag] p1      = (" << winP1.x()        << ", " << winP1.y()        << ", " << winP1.z()        << ")\n"
            << "[sweep-diag] p2      = (" << winP2.x()        << ", " << winP2.y()        << ", " << winP2.z()        << ")\n"
            << "[sweep-diag] p3      = (" << winP3.x()        << ", " << winP3.y()        << ", " << winP3.z()        << ")\n"
            << "[sweep-diag] end     = (" << cornerMirror.x() << ", " << cornerMirror.y() << ", " << cornerMirror.z() << ")\n"
            << "[sweep-diag] T       = " << tr.second << "\n";
  std::cout << std::setprecision(6);   // restore for subsequent output

  // Write the winning trajectory's K+1 waypoints to best_diag_spline.csv (FIX 1).
  {
    std::ofstream splineFile("best_diag_spline.csv");
    splineFile << "k,t,x,y,z\n";
    for (size_t wi = 0; wi < tr.first.size(); ++wi) {
      const TrajectoryPoint& tp = tr.first[wi];
      splineFile << std::setprecision(9) << wi << ',' << tp.time << ','
                 << tp.position.x() << ',' << tp.position.y() << ',' << tp.position.z() << '\n';
    }
  }

  ParsedArgs pa = base;
  pa.cfg = c;
  pa.cfg.customTraj = nullptr;   // printCsv only reads scalar params
  pa.foldName = "diag";
  pa.gripName = "corner";
  printCsv(pa, runFold(c));
  return 0;
}

// ===========================================================================
// --sweep : for the SYMMETRIC fold this is EXECUTION CALIBRATION — LHS(128)
// seeding + diagonal-Gaussian CEM over the 5 execution scalars
// (T, gain, maxForce, eps, delta) (§10); gain/maxForce in log10-space. Cost:
//   J = foldErrMean + 0.5*interiorStrainPeak
//       + 1000*(interpenetrations > 0 || !finite || interiorStrainPeak > 0.15)
// For the DIAGONAL fold, --sweep dispatches to the spline + full-covariance
// CEM path (CHANGE 8), the genuine trajectory optimization. Deterministic.
// ===========================================================================
static int runSweep(const ParsedArgs& base) {
  if (base.cfg.fold == FoldType::Diag) return runSweepDiagonal(base);   // CHANGE 8c

  struct Dim { const char* name; double lo, hi; bool logScale; };
  const std::array<Dim, 5> dims = {{ { "T", 2.0, 10.0, false },
                                     { "gain", 500.0, 8000.0, true },
                                     { "maxForce", 100.0, 800.0, true },
                                     { "eps", 0.01, 0.05, false },
                                     { "delta", 0.0, 0.06, false } }};
  using X = std::array<double, 5>;   // z in [0,1]^5 (log dims mapped through log10)

  auto toValue = [&](int d, double z) {
    const Dim& D = dims[size_t(d)];
    if (D.logScale) return std::pow(10.0, std::log10(D.lo) + z * (std::log10(D.hi) - std::log10(D.lo)));
    return D.lo + z * (D.hi - D.lo);
  };
  auto toConfig = [&](const X& z) {
    RunConfig c = base.cfg;
    c.trace = false;
    c.P.T = toValue(0, z[0]);
    c.P.gain = toValue(1, z[1]);
    c.P.maxForce = toValue(2, z[2]);
    c.P.eps = toValue(3, z[3]);
    c.P.delta = toValue(4, z[4]);
    return c;
  };
  auto cost = [&](const X& z, FoldMetrics& m) {
    m = runFold(toConfig(z));
    double J = m.foldErrMean + 0.5 * m.interiorStrainPeak;
    if (m.interpenetrations > 0 || !m.finite || m.interiorStrainPeak > 0.15) J += 1000.0;
    return J;
  };
  auto show = [&](const char* tag, const X& z, double J) {
    std::cout << tag << " J=" << std::setprecision(5) << J;
    for (int d = 0; d < 5; ++d) std::cout << ' ' << dims[size_t(d)].name << '=' << toValue(d, z[d]);
    std::cout << std::endl;
  };

  std::mt19937 rng(42);
  std::uniform_real_distribution<double> uni(0.0, 1.0);

  // Stage 1: Latin hypercube, 128 samples — one stratum per sample per dim.
  const int NS = 128;
  std::array<std::vector<int>, 5> perms;
  for (auto& p : perms) {
    p.resize(size_t(NS));
    for (int k = 0; k < NS; ++k) p[size_t(k)] = k;
    std::shuffle(p.begin(), p.end(), rng);
  }
  std::vector<std::pair<double, X>> pool;
  for (int k = 0; k < NS; ++k) {
    X z;
    for (int d = 0; d < 5; ++d) z[size_t(d)] = (perms[size_t(d)][size_t(k)] + uni(rng)) / NS;
    FoldMetrics m;
    pool.emplace_back(cost(z, m), z);
  }
  std::sort(pool.begin(), pool.end(),
            [](const auto& a, const auto& b) { return a.first < b.first; });
  pool.resize(16);                                     // keep best 16
  show("[sweep] LHS best:", pool[0].second, pool[0].first);

  // Stage 2: CEM — population 32, elite 8, <= 12 iterations, diagonal
  // Gaussian in z-space with a sigma floor, clipped to bounds.
  X mean{}, sigma{};
  for (int d = 0; d < 5; ++d) {
    double mu = 0.0;
    for (const auto& e : pool) mu += e.second[size_t(d)];
    mu /= double(pool.size());
    double var = 0.0;
    for (const auto& e : pool) var += (e.second[size_t(d)] - mu) * (e.second[size_t(d)] - mu);
    mean[size_t(d)] = mu;
    sigma[size_t(d)] = std::max(std::sqrt(var / double(pool.size())), 0.05);
  }
  std::pair<double, X> best = pool[0];
  std::normal_distribution<double> gauss(0.0, 1.0);
  for (int iter = 1; iter <= 12; ++iter) {
    std::vector<std::pair<double, X>> gen;
    for (int k = 0; k < 32; ++k) {
      X z;
      for (int d = 0; d < 5; ++d)
        z[size_t(d)] = std::clamp(mean[size_t(d)] + sigma[size_t(d)] * gauss(rng), 0.0, 1.0);
      FoldMetrics m;
      gen.emplace_back(cost(z, m), z);
    }
    std::sort(gen.begin(), gen.end(),
              [](const auto& a, const auto& b) { return a.first < b.first; });
    if (gen[0].first < best.first) best = gen[0];
    for (int d = 0; d < 5; ++d) {                     // refit on the 8 elites
      double mu = 0.0;
      for (int e = 0; e < 8; ++e) mu += gen[size_t(e)].second[size_t(d)];
      mu /= 8.0;
      double var = 0.0;
      for (int e = 0; e < 8; ++e)
        var += (gen[size_t(e)].second[size_t(d)] - mu) * (gen[size_t(e)].second[size_t(d)] - mu);
      mean[size_t(d)] = mu;
      sigma[size_t(d)] = std::max(std::sqrt(var / 8.0), 0.01);   // sigma floor
    }
    show(("[sweep] iter " + std::to_string(iter) + " best:").c_str(), best.second, best.first);
  }

  std::cout << "[sweep] FINAL" << std::endl;
  ParsedArgs pa = base;
  pa.cfg = toConfig(best.second);
  pa.cfg.trace = base.cfg.trace;
  printCsv(pa, runFold(pa.cfg));
  return 0;
}

// ===========================================================================
// GUI mode — modeled on optimization.cpp's working idioms.
// ===========================================================================
namespace gui {

struct State {
  int N = 10;
  FoldType fold = FoldType::SymZ;
  FoldParams P;

  ClothMesh cloth;
  ClothController controller;
  std::unique_ptr<SemiImplicitEulerIntegrator> integrator;
  polyscope::SurfaceMesh* psMesh = nullptr;

  bool folding = false;
  double simTime = 0.0;
  bool useGfold = false;   // CHANGE 10: drive the OneCorner diagonal with the g-fold path

  bool haveSpec = false;
  FoldSpec spec;
  std::vector<SpringRef> springs;
  double liveFoldErr = 0.0, liveStrain = 0.0;
};

// Reposition to the flat start pose at y = h; unpin everything (E10: no
// unpin API — write pinned[i] = false directly); clear all controls. The
// controller CLOCK is deliberately NOT reset: trajectories are relative to
// their own start time (audited setTrajectory behavior).
static void reset(State& st) {
  for (size_t v = 0; v < st.cloth.getVertexCount(); ++v) {
    const int i = static_cast<int>(v) / st.N;
    const int j = static_cast<int>(v) % st.N;
    st.cloth.setVertexPosition(v, Eigen::Vector3d(j * kSpacing, st.P.h, i * kSpacing));
    st.cloth.velocities[v] = Eigen::Vector3d::Zero();
    st.cloth.pinned[v] = false;
  }
  st.controller.removeAllControls();
  st.folding = false;
  st.haveSpec = false;
  st.simTime = 0.0;
  st.liveFoldErr = st.liveStrain = 0.0;
}

// Fresh ClothMesh at the current N; grid repositioned to y = h (createGrid
// hardcodes y = 1.0); properties + gravity per E1/E2b; fresh controller
// state; re-register the Polyscope mesh under the same name.
static void rebuild(State& st) {
  st.cloth = ClothMesh();
  st.cloth.createGrid(st.N, st.N, kSpacing);
  reset(st);
  st.cloth.properties.stiffness = st.P.stiffness;
  st.cloth.properties.bendingStiffness = st.P.bending;
  st.cloth.properties.damping = st.P.damping;
  st.cloth.properties.gravity = resolveGravityScale(st.P, st.N) * Eigen::Vector3d(0.0, -9.81, 0.0);
  st.psMesh = polyscope::registerSurfaceMesh("Cloth", st.cloth.getVertexMatrix(),
                                             st.cloth.getTriangleMatrix());
  st.psMesh->setSurfaceColor({ 0.2, 0.8, 0.4 });
  st.psMesh->setEdgeWidth(1.0);
}

// Reset first; compute spec; pin S; assert grips disjoint from pins (E10);
// generate arcs and hand them to the controller via the E6 choke point.
static void startFold(State& st) {
  reset(st);
  FoldingOptimizer opt(st.N, st.N, kSpacing);
  RunConfig cfg;
  cfg.fold = st.fold;
  cfg.N = st.N;
  cfg.P = st.P;
  const FoldSpec spec = makeSpec(cfg, opt);
  if (spec.moving.empty()) {
    std::cerr << "[Folding] empty fold spec; not starting." << std::endl;
    return;
  }
  st.springs = classifySprings(st.cloth, spec);

  for (size_t v : spec.pinnedSet) st.cloth.pinVertex(static_cast<int>(v));

  const std::vector<size_t> grips = opt.gripSet(spec, st.P.grip);
  for (size_t g : grips) {
    if (st.cloth.pinned[g]) {   // E10: controlled ∩ pinned must be empty
      std::cerr << "[Folding] FATAL: grip vertex " << g << " is pinned (E10)." << std::endl;
      std::abort();
    }
  }

  std::vector<TrajectoryPoint> traj;
  FoldParams params = st.P;
  if (st.useGfold && st.fold == FoldType::Diag)   // CHANGE 10: g-fold gravity path
    traj = opt.gFoldTrajectory(st.cloth, Corner::FarRight, Corner::NearLeft, params);
  else if (st.fold == FoldType::Diag)
    traj = opt.optimizeDiagonalFold(st.cloth, Corner::FarRight, Corner::NearLeft, st.P.T, params);
  else
    traj = opt.optimizeSymmetricFold(st.cloth, st.fold == FoldType::SymZ ? FoldDirection::AlongX
                                                                         : FoldDirection::AlongZ,
                                     st.P.T, params);
  opt.applyFoldingTrajectory(st.controller, traj, st.P.gain, st.P.maxForce);

  st.spec = spec;
  st.haveSpec = true;
  st.folding = true;
  st.simTime = 0.0;
  std::cout << "[Folding] fold started: " << traj.size() << " waypoints across "
            << grips.size() << " gripped vertices." << std::endl;
}

static int run() {
  polyscope::options::autocenterStructures = false;
  polyscope::options::autoscaleStructures = false;
  polyscope::view::windowWidth = 1600;
  polyscope::view::windowHeight = 1000;
  polyscope::init();
  polyscope::view::setUpDir(polyscope::UpDir::YUp);
  polyscope::view::setFrontDir(polyscope::FrontDir::ZFront);

  State st;
  st.integrator = std::make_unique<SemiImplicitEulerIntegrator>();
  rebuild(st);

  // Visual-only ground quad at y = 0 (the integrator clamps there).
  {
    std::vector<Eigen::Vector3d> gv = { { -0.4, 0.0, -0.4 }, { 1.4, 0.0, -0.4 },
                                        { 1.4, 0.0, 1.4 },  { -0.4, 0.0, 1.4 } };
    std::vector<std::array<int, 3>> gt = { { 0, 1, 2 }, { 0, 2, 3 } };
    auto* ground = polyscope::registerSurfaceMesh("Ground", gv, gt);
    ground->setSurfaceColor({ 0.7, 0.7, 0.7 });
    ground->setMaterial("flat");
  }
  polyscope::view::lookAt({ 1.8, 1.2, 1.8 }, { 0.45, 0.15, 0.45 }, { 0.0, 1.0, 0.0 });

  polyscope::state::userCallback = [&]() {
    // Live-apply the physics sliders + gravity per E1/E2b every frame so
    // toggling realGravity mid-run behaves.
    st.cloth.properties.stiffness = st.P.stiffness;
    st.cloth.properties.bendingStiffness = st.P.bending;
    st.cloth.properties.damping = st.P.damping;
    st.cloth.properties.gravity = resolveGravityScale(st.P, st.N) * Eigen::Vector3d(0.0, -9.81, 0.0);

    // Demo loop order: controls -> physics -> clamps -> render update. The
    // integrator is subdivided into P.substeps sub-steps (CHANGE 4b) so the
    // GUI matches the headless rollout at high stiffness / restored gravity.
    if (st.folding) {
      st.controller.applyControls(st.cloth, st.P.dt);
      st.simTime += st.P.dt;
    }
    const int guiSubsteps = std::max(1, st.P.substeps);
    for (int sub = 0; sub < guiSubsteps; ++sub) st.integrator->step(st.cloth, st.P.dt / guiSubsteps);

    // E8: the Optimization demo's soft ground clamp, replicated verbatim.
    for (size_t i = 0; i < st.cloth.getVertexCount(); ++i) {
      Eigen::Vector3d p = st.cloth.getVertex(i).position;
      if (p.y() < 0.01) {
        p.y() = 0.01;
        st.cloth.setVertexPosition(i, p);
        st.cloth.velocities[i].y() *= 0.3;
      }
    }

    // Virtual layer clamp (§10): stand-in for missing layer contact,
    // generalizing the engine's own ground-clamp idiom to the pinned half's
    // plane — moving vertices past the crease may not sink below h + eps.
    if (st.P.layerClamp && st.haveSpec) {
      for (size_t v : st.spec.moving) {
        Eigen::Vector3d p = st.cloth.getVertex(v).position;
        const double rNow = (p - st.spec.creasePoint).dot(st.spec.creaseNormal);
        if (rNow < -kTol && p.y() < st.P.h + st.P.eps) {
          p.y() = st.P.h + st.P.eps;
          st.cloth.setVertexPosition(v, p);
          st.cloth.velocities[v].y() *= 0.3;
        }
      }
    }

    // Live metrics: current M1 fold error and current interior strain.
    if (st.haveSpec) {
      const double L = (st.N - 1) * kSpacing;
      double sum = 0.0;
      for (size_t k = 0; k < st.spec.moving.size(); ++k)
        sum += (st.cloth.getVertex(st.spec.moving[k]).position - st.spec.mirrorTargets[k]).norm() / L;
      st.liveFoldErr = sum / double(st.spec.moving.size());
      double peak = 0.0;
      for (const SpringRef& sp : st.springs) {
        if (sp.crossing) continue;
        const double len = (st.cloth.getVertex(size_t(sp.v0)).position -
                            st.cloth.getVertex(size_t(sp.v1)).position).norm();
        peak = std::max(peak, std::abs(len - sp.rest) / sp.rest);
      }
      st.liveStrain = peak;
    }

    st.psMesh->updateVertexPositions(st.cloth.getVertexMatrix());

    // ------------------------- GUI panel -------------------------
    ImGui::SetNextWindowSize(ImVec2(430, 780), ImGuiCond_FirstUseEver);
    ImGui::SetNextWindowPos(ImVec2(40, 40), ImGuiCond_FirstUseEver);
    if (ImGui::Begin("Folding Optimizer", nullptr, ImGuiWindowFlags_None)) {
      ImGui::Text("=== GRID ===");
      ImGui::SliderInt("N (grid)", &st.N, 5, 20);
      ImGui::SameLine();
      if (ImGui::Button("Rebuild")) rebuild(st);

      ImGui::Separator();
      ImGui::Text("=== FOLD TYPE ===");
      int fold = static_cast<int>(st.fold);
      if (ImGui::RadioButton("Symmetric Z", &fold, 0)) st.P.T = 5.0;
      if (ImGui::RadioButton("Symmetric X", &fold, 1)) st.P.T = 5.0;
      if (ImGui::RadioButton("Diagonal FarRight->NearLeft", &fold, 2)) st.P.T = 6.0;  // §7
      st.fold = static_cast<FoldType>(fold);

      ImGui::Text("=== GRIP ===");
      int grip = static_cast<int>(st.P.grip);
      ImGui::RadioButton("All", &grip, 0); ImGui::SameLine();
      ImGui::RadioButton("FarEdge", &grip, 1); ImGui::SameLine();
      ImGui::RadioButton("TwoCorners", &grip, 2); ImGui::SameLine();
      ImGui::RadioButton("OneCorner", &grip, 3);
      st.P.grip = static_cast<GripMode>(grip);

      ImGui::Separator();
      ImGui::Text("=== PARAMETERS ===");
      ImGui::PushItemWidth(220);
      auto sliderD = [](const char* label, double* v, float lo, float hi, const char* fmt) {
        float f = static_cast<float>(*v);
        if (ImGui::SliderFloat(label, &f, lo, hi, fmt)) *v = double(f);
      };
      sliderD("T (fold duration)", &st.P.T, 2.0f, 10.0f, "%.1f");
      sliderD("gain", &st.P.gain, 500.0f, 8000.0f, "%.0f");
      sliderD("maxForce", &st.P.maxForce, 100.0f, 800.0f, "%.0f");
      sliderD("eps (layer gap)", &st.P.eps, 0.005f, 0.05f, "%.3f");
      sliderD("delta (crease slack)", &st.P.delta, 0.0f, 0.06f, "%.3f");
      // CHANGE 9: elliptical-family + stability knobs. alpha shapes the arc
      // height (1 = rigid arc); substeps and gravityScale are the stability /
      // drape knobs. The GUI reads st.P live into the arc + cloth each frame,
      // so these flow through automatically. Run --grid / --sweep headless,
      // then set the sliders to the printed best (optimizer is NOT in the GUI).
      sliderD("alpha (arc height)", &st.P.alpha, 0.5f, 1.2f, "%.3f");
      sliderD("gravityScale", &st.P.gravityScale, 1.0f, 200.0f, "%.0f");
      ImGui::SliderInt("substeps", &st.P.substeps, 1, 8);
      sliderD("stiffness", &st.P.stiffness, 100.0f, 5000.0f, "%.0f");
      sliderD("damping", &st.P.damping, 0.3f, 0.95f, "%.2f");
      ImGui::PopItemWidth();
      ImGui::Checkbox("realGravity (E2b: g x N^2)", &st.P.realGravity);
      ImGui::SameLine();
      ImGui::Checkbox("layerClamp", &st.P.layerClamp);
      ImGui::Checkbox("g-fold path (OneCorner diagonal)", &st.useGfold);

      ImGui::Separator();
      ImGui::PushStyleColor(ImGuiCol_Button, ImVec4(0.2f, 0.7f, 0.2f, 1.0f));
      if (ImGui::Button("Start Fold", ImVec2(120, 0))) startFold(st);
      ImGui::PopStyleColor();
      ImGui::SameLine();
      if (ImGui::Button("Reset", ImVec2(100, 0))) reset(st);
      ImGui::SameLine();
      ImGui::PushStyleColor(ImGuiCol_Button, ImVec4(0.8f, 0.6f, 0.2f, 1.0f));
      if (ImGui::Button("Drape preset", ImVec2(140, 0))) {
        // §7: real gravity organizes the slack; the layer clamp catches it.
        st.P.realGravity = true;
        st.P.layerClamp = true;
        st.P.stiffness = 2500.0;
        st.P.settle = 2.5;
        st.P.substeps = 8;
      }
      ImGui::PopStyleColor();

      ImGui::Separator();
      ImGui::Text("=== STATUS ===");
      ImGui::TextColored(st.folding ? ImVec4(0.2f, 1.0f, 0.2f, 1.0f)
                                    : ImVec4(1.0f, 0.4f, 0.4f, 1.0f),
                         "Folding: %s", st.folding ? "ACTIVE" : "idle");
      ImGui::Text("Sim time: %.2f s", st.simTime);
      ImGui::Text("Controlled vertices: %zu", st.controller.getControlCount());
      ImGui::Text("alpha: %.3f%s", st.P.alpha, st.useGfold ? "  [g-fold]" : "");
      ImGui::Text("Live fold error D (M1): %.4f of L", st.liveFoldErr);
      ImGui::Text("Live interior strain S (M2): %.4f", st.liveStrain);
    }
    ImGui::End();
  };

  polyscope::show();
  return 0;
}

} // namespace gui

int main(int argc, char** argv) {
  // §4 invariant checks first; stderr keeps --eval's stdout CSV clean.
  if (!selfTest()) return 1;

  bool eval = false, sweep = false, grid = false, gfold = false;
  for (int i = 1; i < argc; ++i) {
    if (std::strcmp(argv[i], "--eval") == 0) eval = true;
    if (std::strcmp(argv[i], "--sweep") == 0) sweep = true;
    if (std::strcmp(argv[i], "--grid") == 0) grid = true;
    if (std::strcmp(argv[i], "--gfold") == 0) gfold = true;
  }

  if (eval || sweep || grid || gfold) {
    ParsedArgs pa;
    if (!parseArgs(argc, argv, pa)) return 2;
    if (sweep) return runSweep(pa);      // symmetric CEM or diagonal spline+CMA (CHANGE 8)
    if (grid) {
      if (pa.cfg.fold == FoldType::Diag) return runGridDiagonal(pa);
      return runGrid(pa);               // exhaustive alpha x T (CHANGE 7)
    }
    return runEval(pa);                  // --eval, or --gfold single run (CHANGE 10)
  }
  return gui::run();
}

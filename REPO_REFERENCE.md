# cloth_opt — Repository Reference

**Provenance:** branch `tiny_opt`, commit `62bf397c0174b47130865478c392ea3d962c4c40` ("fix gitignore bug", 2025-06-14). Working tree clean at time of writing.

**Method note:** Everything below was determined by reading the source (static analysis). The project was **not built or run** on this machine (no `build/` directory exists; the build requires system OpenVDB + LibXml2, see §8). Statements that are deterministic consequences of the code (formulas, constants) are stated as fact; anything requiring execution to confirm is labeled **[inference]** or listed in §11.

## Table of Contents

1. [Overview and purpose](#1-overview-and-purpose)
2. [Directory and file map](#2-directory-and-file-map)
3. [Data model](#3-data-model)
4. [Coordinate system and world layout](#4-coordinate-system-and-world-layout)
5. [Physics and integrator](#5-physics-and-integrator)
6. [Control system](#6-control-system) ← most important
7. [The three demos](#7-the-three-demos)
8. [Build system](#8-build-system)
9. [Public API quick reference](#9-public-api-quick-reference)
10. [Notes for the folding task](#10-notes-for-the-folding-task)
11. [Open questions and runtime unknowns](#11-open-questions-and-runtime-unknowns)

---

## 1. Overview and purpose

ClothOpt is a small C++17 mass-spring cloth simulator with an interactive per-vertex control system, rendered with Polyscope. All project code lives in namespace `ClothOpt`. The architecture:

- **`SurfaceMesh`** (`include/mesh.h`, `src/mesh.cpp`) — plain triangle mesh: vertices (position + normal), triangles, Eigen matrix export.
- **`ClothMesh`** (`include/cloth.h`, `src/cloth.cpp`) — extends `SurfaceMesh` with per-vertex simulation state (velocity, force, pinned flag), spring constraints (structural + bending; **no shear**), collision spheres, and `ClothProperties`.
- **`SemiImplicitEulerIntegrator`** (`include/integrator.h`, `src/integrator.cpp`) — one concrete integrator: gravity + spring forces → semi-implicit Euler with multiplicative velocity damping → sphere/ground collision resolution → 1 iteration of position-based distance-constraint projection.
- **`ClothController`** (`include/controller.h`, `src/controller.cpp`) — per-vertex P/D/force controllers plus trajectory following and circular/sinusoidal motion patterns. **Applies control by directly adding Δv to `mesh.velocities`**, not via the force array.
- **Three demos** (top-level `.cpp` files → CMake targets): `Visualization` (static wave-deformed grid, does not use ClothOpt classes at all), `Simulation` (8×8 cloth dropped onto a sphere, no controller), `Optimization` (on this branch: interactive 2×2 cloth with one GUI-controlled vertex; on `master` it is a 10×10 edge-pulling demo).

Typical frame loop (from `optimization.cpp:126-142`): `controller.applyControls(cloth, dt)` → `integrator->step(cloth, dt)` → ad-hoc ground clamp → update Polyscope.

## 2. Directory and file map

### Project code

| File | Responsibility |
|---|---|
| `include/mesh.h`, `src/mesh.cpp` | `Vertex`, `Triangle`, `SurfaceMesh` (add/clear/get/set, `computeNormals`, `getVertexMatrix`/`getTriangleMatrix`), `MeshUtils` (triangle normal/area). |
| `include/cloth.h`, `src/cloth.cpp` | `Edge`, `CollisionSphere`, `ClothProperties`, `ClothMesh` (`createGrid`, `buildConstraints`, `pinVertex`, `pinCorners`, `addSphere`, `clearSpheres`, `getGridIndex`). |
| `include/integrator.h`, `src/integrator.cpp` | `DebugInfo`, abstract `Integrator`, `SemiImplicitEulerIntegrator` (`step`, gravity/springs, PBD constraint projection, sphere+ground collisions, `pointToTriangleDistance`, debug printers). Contains a **dead** `updatePositions` function (never called; see §5). |
| `include/controller.h`, `src/controller.cpp` | `ControlTarget`, `ClothController` (position/velocity/force controls, trajectories, motion patterns, wind). |
| `visual.cpp` | Target `Visualization`. Pure-Polyscope static demo; builds its own 20×20 grid + sine-wave deformation. Uses **none** of the ClothOpt classes. |
| `simulation.cpp` | Target `Simulation`. 8×8 ClothMesh falls onto a collision sphere; two pinned vertices; performance instrumentation; no controller. |
| `optimization.cpp` | Target `Optimization`. **This branch:** 2×2 cloth, vertex 0 pinned, GUI drives one vertex via `addPositionControl`. |
| `CMakeLists.txt` | Build config; three executables, each compiled with all of `src/*.cpp` (no library target). |
| `README.md` | Setup + docs. Several claims are stale on this branch — flagged throughout as **[README≠code]**. |
| `image/` | Screenshots for README. `.vscode/` — a generic single-file g++ task (not the real build). |

### `deps/` survey (vendored; internals not audited)

| Dep | What it is | Actually used by project code? |
|---|---|---|
| `deps/polyscope` | 3D viewer + bundled ImGui | **Yes** — all three demos render with it; `optimization.cpp` uses ImGui directly for its GUI. Built via `add_subdirectory` (`CMakeLists.txt:123`). |
| `deps/eigen` | Eigen (vendored copy) | **Yes (headers)** — all math is Eigen. Note the build both links system `Eigen3::Eigen` (`CMakeLists.txt:65,169`) *and* adds `deps/eigen` to include paths (`CMakeLists.txt:160`); which copy's headers win depends on include order. |
| `deps/alembic` | Alembic scene-interchange lib | **No** — linked (`CMakeLists.txt:171`) but no project source includes it. Presumably intended for future export. |
| `deps/Imath` | Math lib required by Alembic | **No** — built (`CMakeLists.txt:126`) only as Alembic's dependency. |
| `deps/openvdb` | OpenVDB (vendored copy) | **No, and not even built** — CMake instead does `find_package(OpenVDB REQUIRED)` against the **system** (`CMakeLists.txt:118-119`) and links it; no project source uses it. The vendored copy is dead weight. |
| `deps/libigl` | libigl geometry headers | **No** — include path only (`CMakeLists.txt:155`); never included. |
| `deps/json` | nlohmann/json | **No** — include path only (`CMakeLists.txt:161`); never included. |
| `deps/fbx` | Referenced at `CMakeLists.txt:156` | **Directory does not exist.** Harmless stale include path. |

## 3. Data model

### Grid creation — `ClothMesh::createGrid(int width, int height, double spacing = 0.1)` (`src/cloth.cpp:9-44`)

- `width` = number of **columns** (index `j`), `height` = number of **rows** (index `i`).
- Vertex `(i, j)` is created at world position **`(x, y, z) = (j*spacing, 1.0, i*spacing)`** (`src/cloth.cpp:18`). The initial height `y = 1.0` is hardcoded; every demo immediately overwrites `y` afterwards.
- Creation order: row-major, `i` outer, `j` inner ⇒ **linear index = `i * width + j`**.
- Public helper: **`int getGridIndex(int i, int j) const { return i * gridWidth_ + j; }`** (`include/cloth.h:66`). Note argument order: **row first, column second**.
- Triangles per cell: `(topLeft, bottomLeft, topRight)` and `(topRight, bottomLeft, bottomRight)` (`src/cloth.cpp:31-32`), which yields upward (+y) normals for the initial flat grid.
- After vertices/triangles: `velocities`, `forces`, `pinned` are resized to vertex count and zero/false-initialized; `buildConstraints()`; `computeNormals()` (`src/cloth.cpp:37-43`).

### Per-vertex state

| Field | Where | Notes |
|---|---|---|
| `position`, `normal` | `Vertex` in `SurfaceMesh::vertices_` (`include/mesh.h:12-18`) | Read via `getVertex(i).position`; write via `setVertexPosition(i, p)` (silently ignores out-of-range index, `src/mesh.cpp:21-25`). |
| `velocities[i]` | `ClothMesh` public `std::vector<Eigen::Vector3d>` | Directly writable (the controller does). |
| `forces[i]` | same | Cleared at the start of every `step`. The controller does **not** use it. |
| `pinned[i]` | `std::vector<bool>` | Public. `pinVertex` sets it true + zeroes velocity (`src/cloth.cpp:97-102`). **There is no unpin API** — write `mesh.pinned[i] = false;` directly. |

### Constraints (`ClothMesh::buildConstraints`, `src/cloth.cpp:46-87`)

- `Edge {int v0, v1; double restLength;}` (`include/cloth.h:11-15`).
- **`distanceConstraints` (structural):** every horizontal neighbor `(i,j)-(i,j+1)` and vertical neighbor `(i,j)-(i+1,j)`, rest length = `spacing`. Count = `H*(W-1) + W*(H-1)`.
- **`bendingConstraints`:** skip-one springs `(i,j)-(i,j+2)` and `(i,j)-(i+2,j)`, rest length = `2*spacing`.
- **There are NO shear (diagonal) constraints of any kind.** Triangles are render-only geometry. In-plane shearing of a grid cell is resisted by nothing (the bending springs are collinear with the structural directions). **[README≠code]** — irrelevant to the README but a key physical fact.
- `collisionSpheres`: list of `{center, radius}` added via `addSphere` (`src/cloth.cpp:89-95`).

## 4. Coordinate system and world layout

- **Y is up.** Ground plane is `y = 0`, enforced inside the integrator (`src/integrator.cpp:232-236`: any non-pinned vertex with `position.y() < 0` is clamped to `y = 0`). Both simulating demos also call `polyscope::view::setUpDir(YUp)`.
- **Grid rows (`i`) run along +Z; grid columns (`j`) run along +X; the sheet is horizontal (constant y).**
- **Gravity mismatch:** `ClothProperties` **default** is `gravity = (0, 0, -9.81)` (`include/cloth.h:31`) — i.e. along **−Z**, which is *sideways* relative to the y=0 ground. Both simulating demos override it to `(0, -9.81, 0)` (`simulation.cpp:50`, `optimization.cpp:56`). **If you write a new demo and forget to set gravity, the cloth accelerates horizontally.**
- Example — `createGrid(10, 10, 0.1)` (as on master), before any demo repositioning (y = 1.0 from `createGrid`):

| Corner | (i, j) | linear index | world position |
|---|---|---|---|
| origin corner | (0,0) | 0 | (0.0, 1.0, 0.0) |
| +x corner | (0,9) | 9 | (0.9, 1.0, 0.0) |
| +z corner | (9,0) | 90 | (0.0, 1.0, 0.9) |
| far corner | (9,9) | 99 | (0.9, 1.0, 0.9) |

- `pinCorners()` pins exactly these four indices: `0`, `W−1`, `(H−1)*W`, `H*W−1` (`src/cloth.cpp:104-111`). No demo currently calls it.

## 5. Physics and integrator

`SemiImplicitEulerIntegrator::step(ClothMesh&, double dt)` (`src/integrator.cpp:20-54`) does, in order:

1. **Zero `mesh.forces`.**
2. **Gravity** (`applyGravity`, `:56-82`): per non-pinned vertex adds `(properties.mass / vertexCount) * gravity` to `forces[i]`.
3. **Spring forces** (`applySpringForces`, `:84-135`): for each structural edge, `F = stiffness * (len − rest) * dir` applied ± to the two endpoints; same for bending edges with `bendingStiffness`. No spring damping term. Zero-length guard at 1e-8.
4. **Velocity + position update** (inline in `step`, `:31-41`), skipping pinned vertices:
   - `v ← v * damping + (F / properties.mass) * dt`  (damping multiplies **old** velocity only)
   - `x ← x + v * dt`
5. **Collisions** (`handleCollisions`, `:214-399`) — position projection only, three TBB-parallel passes:
   - **Vertex:** ground clamp `y<0 → y=0`; vertex-in-sphere pushed to sphere surface.
   - **Edge–sphere:** closest point on each structural edge vs each sphere; both endpoints pushed out by half the penetration.
   - **Face–sphere:** `pointToTriangleDistance` (`:402-490`); all 3 vertices pushed out by penetration/3.
   - **No velocity change on contact** (no restitution, no friction, velocities keep pointing into the ground/sphere). **`properties.friction` is never read anywhere in the codebase** — dead parameter.
   - **No self-collision and no cloth–cloth collision of any kind.** Only sphere and ground-plane handling exist.
6. **Constraint projection** (`satisfyConstraints`, `:186-212`): **exactly 1 iteration** (`const int iterations = 1;`, `:187`) of PBD-style projection over **structural constraints only** (bending constraints are never projected). Each endpoint moves half the correction unless pinned. Note this runs **after** collisions, so projection can push vertices back below ground / into the sphere within a step.

### Gotchas / quirks (all static facts)

- **Effective gravity is divided by vertex count.** Gravity force per vertex is `(mass/N)*g` (`:57`), but the velocity update divides by **total** `properties.mass` (`:35`), so gravitational acceleration = `g/N`. With defaults: 2×2 grid → 2.45 m/s²; 8×8 → 0.153 m/s²; 10×10 → 0.098 m/s². Spring accelerations are *not* scaled by N. Consequence: cloth falls much slower than real gravity and the effect worsens with resolution. **[inference about intent: this looks like a bug, mass bookkeeping is inconsistent.]**
- **`updatePositions` (`:137-184`) is dead code** — the physically consistent version (divides by per-vertex mass, damps after adding acceleration) that `step` never calls.
- `handleCollisions` uses `tbb::parallel_for`/`tbb::combinable` **unconditionally** (`:215-217` etc.), even though the `#include <tbb/...>` block is guarded by `#ifdef HAVE_TBB` (`:6-16`, whose comment claims "not used"). **Without TBB the file does not compile** — TBB is a hard dependency in practice (see §8).
- **[inference]** The parallel collision passes write shared vertices from multiple threads (adjacent edges/faces share vertices; `setVertexPosition` is unsynchronized) — a data race in principle; practical impact unknown without running.
- Damping is per-step multiplicative (default 0.99 ⇒ ~0.6× velocity retained per 50 steps), so behavior depends on dt choice.

### Tunable parameters — `ClothProperties` defaults (`include/cloth.h:25-32`)

| Parameter | Default | Used where |
|---|---|---|
| `mass` | 1.0 | gravity (÷N) and velocity update (total) |
| `stiffness` | 1000.0 | structural spring force |
| `bendingStiffness` | 100.0 | bending spring force |
| `damping` | 0.99 | per-step velocity multiplier |
| `friction` | 0.8 | **never used** |
| `gravity` | (0, 0, −9.81) | **wrong axis vs. ground; every demo overrides to (0, −9.81, 0)** |

## 6. Control system

Header `include/controller.h`, implementation `src/controller.cpp`. One `ControlTarget` per controlled vertex, stored in `unordered_map<size_t, ControlTarget>`. A single target can simultaneously have position, velocity, and force control active — **they share one `maxForce` field** (see gotchas).

### `ControlTarget` defaults (`include/controller.h:24-34`)

`positionGain = 1000`, `velocityGain = 100`, `maxForce = 100`.

### Public methods (exact signatures, `include/controller.h`)

```cpp
void addPositionControl(size_t vertexIndex, const Eigen::Vector3d& targetPosition,
                        double gain = 1000.0, double maxForce = 100.0);
void addVelocityControl(size_t vertexIndex, const Eigen::Vector3d& targetVelocity,
                        double gain = 100.0, double maxForce = 50.0);
void addForceControl(size_t vertexIndex, const Eigen::Vector3d& force);
void updatePositionTarget(size_t vertexIndex, const Eigen::Vector3d& newTarget);   // no-op if vertex not already controlled
void updateVelocityTarget(size_t vertexIndex, const Eigen::Vector3d& newVelocity); // "
void updateForceTarget(size_t vertexIndex, const Eigen::Vector3d& newForce);       // "
void removeControl(size_t vertexIndex);        // also erases trajectory + motion pattern
void removeAllControls();
void applyControls(ClothMesh& mesh, double dt);
void setTrajectory(size_t vertexIndex, const std::vector<Eigen::Vector3d>& positions,
                   const std::vector<double>& times, bool loop = false);
void updateTrajectories(double currentTime);   // called internally by applyControls
bool isControlled(size_t vertexIndex) const;
Eigen::Vector3d getControlForce(size_t vertexIndex) const;  // returns ONLY externalForce, NOT the PD force
size_t getControlCount() const;
std::vector<size_t> getControlledVertices() const;
void enableDebug(bool enable);
void addCircularMotion(size_t vertexIndex, const Eigen::Vector3d& center,
                       double radius, double frequency,
                       const Eigen::Vector3d& axis = Eigen::Vector3d(0,1,0));
void addSinusoidalMotion(size_t vertexIndex, const Eigen::Vector3d& center,
                         const Eigen::Vector3d& amplitude, double frequency);
void addWindForce(const std::vector<size_t>& vertexIndices,
                  const Eigen::Vector3d& windDirection, double strength, double turbulence = 0.0);
```

### What `applyControls(mesh, dt)` does each call (`src/controller.cpp:104-155`)

1. `currentTime_ += dt` — **the controller has its own clock that only advances when `applyControls` is called.**
2. Motion patterns and trajectories update each affected target's `targetPosition`.
3. For every `ControlTarget` (skipping out-of-range vertex indices):
   - `F = 0`
   - position control active: `F += positionGain * (targetPosition − currentPos)` — **pure P, no derivative term** (`:272-276`)
   - velocity control active: `F += velocityGain * (targetVelocity − currentVel)` (`:278-282`)
   - force control active: `F += externalForce`
   - **clamp:** `if (‖F‖ > maxForce) F = F.normalized() * maxForce` (`:139-141`)
   - **apply:** `mesh.velocities[i] += F * dt` (`:144-148`) — comment says "Assume unit mass". So control **writes velocity directly**; it never touches `mesh.forces`, and it ignores `properties.mass`.

### Force law summary

- Saturation position error = `maxForce / positionGain`. **Defaults: 100/1000 = 0.1 m** — beyond 0.1 m of error the pull force is constant at 100 N-equivalent. In the Optimization demo (gain=slider default 800, maxForce hardcoded 100): saturates at 0.125 m.
- Max velocity change per `applyControls` call = `maxForce * dt` (e.g. 100 × 0.01 = 1.0 m/s per frame in the Optimization demo).
- The Δv is applied **before** the integrator step in the demo loop, so it is then multiplied by `damping` and fights gravity/springs/constraint projection inside `step`.

### `setTrajectory` semantics (`src/controller.cpp:159-187`, interpolation `:284-303`)

- Requires `positions.size() == times.size()` and non-empty, else prints an error and does nothing.
- `startTime` is captured as the controller's `currentTime_` at call time; on every `applyControls`, target = `interpolateTrajectory(traj, currentTime_ − startTime)`. Times are therefore **relative controller-seconds from the setTrajectory call** (and only advance while `applyControls` is being called with its dt).
- **Linear interpolation** between consecutive waypoints.
- `loop = true`: time wraps via `fmod(time, times.back())`.
- After the final time (non-loop), the target holds at `positions.back()` forever.
- **Gotcha 1 (silent gain reset):** `setTrajectory` internally calls `addPositionControl(vertexIndex, positions[0])` with **default arguments** (`:173`), which **overwrites any previously configured `positionGain` back to 1000 and `maxForce` back to 100**. Workaround: call `addPositionControl(v, anything, yourGain, yourMaxForce)` **after** `setTrajectory` — the trajectory only rewrites `targetPosition` each frame, so gains set afterwards stick.
- **Gotcha 2 (pre-start behavior):** if `relativeTime < times[0]`, the segment search falls through and returns `positions.back()` — the **last** waypoint, not the first (`:294-302`). Always make `times[0] = 0.0`. (Same issue recurs after each loop wrap if `times[0] > 0`.)
- **Gotcha 3:** times are assumed strictly increasing; duplicate consecutive times divide by zero **[inference: produces NaN target]**.

### Motion patterns (`src/controller.cpp:213-248`, `:305-331`)

- Circular: target orbits `center` with `radius` at `frequency` Hz around `axis` (default +Y). Basis built as `v = axis × (1,0,0)` (fallback `axis × (0,1,0)`); **latent bug:** if `axis` is parallel to the X axis, the first cross product is the zero vector and `.normalized()` yields NaN, and the `v.norm() < 0.1` fallback check does not catch NaN **[inference: NaN target positions]**. Default axis is safe.
- Sinusoidal: `target = center + amplitude * sin(2π f t)` (component-wise via vector amplitude).
- Both simply rewrite `targetPosition`; also registered via `addPositionControl(...)` with default gain/maxForce (same reset gotcha as trajectories).
- `addWindForce`: registers a **constant** `addForceControl` per vertex; turbulence is a single random sample at call time, not per-step noise.

### Controlled vs pinned vertices

- The integrator **skips pinned vertices entirely** (`src/integrator.cpp:32`) — no forces, no velocity/position update, no collisions, no constraint projection movement.
- The controller does **not** check `pinned`. Controlling a pinned vertex accumulates velocity every frame (`v += F*dt`) that is never applied or damped: **position never moves and the stored velocity grows without bound.** Unpin first (`mesh.pinned[i] = false;` — no API for it), or never control pinned vertices.
- A controlled (non-pinned) vertex is still fully subject to gravity, springs, collisions, and PBD projection; control is just an extra Δv each frame.

## 7. The three demos

### `Visualization` (`visual.cpp`)

- Does **not** use ClothOpt classes; builds its own 20×20 grid (spacing 0.1, y=0) with a static sine-wave deformation (amplitude 0.2, frequency 3.0). No physics, no loop callback, no GUI beyond Polyscope defaults. Purely a rendering sanity check.

### `Simulation` (`simulation.cpp`)

- `createGrid(8, 8, 0.125)`; all vertices then moved to `y = 0.8` (`:35-39`).
- Pins indices `getGridIndex(0,0)`=0 and `getGridIndex(0,7)`=7 (the two corners of the z=0 edge).
- Properties: stiffness 2000, bendingStiffness 200, damping 0.98, friction 0.7 (unused), gravity **(0, −9.81, 0)** (`:46-50`).
- One collision sphere at `(0.4375, 0.5, 0.4375)`, radius 0.12 (cloth center).
- Loop: **8 integrator substeps per rendered frame at dt = 0.003** (`:108-117`); mesh render update only every 2nd frame; console perf stats every 60 frames. No controller, no ad-hoc ground clamp (relies on the integrator's y<0 clamp). Visual-only ground quad at y=0.
- GUI: none (no user callback widgets).

### `Optimization` (`optimization.cpp`) — this branch (`tiny_opt`)

- **`createGrid(2, 2, 0.2)` — a 2×2, 4-vertex cloth** (`:29`). Vertices then placed at `y = 0.05` (`:33-39`). Layout (row i → z, col j → x):
  - v0=(0, .05, 0), v1=(0.2, .05, 0), v2=(0, .05, 0.2), v3=(0.2, .05, 0.2). Console art: `2---3 / 0---1`.
- Properties: stiffness 500, damping 0.8, gravity (0, −9.81, 0). Bending stiffness left at default 100 (only one skip-one pair per direction doesn't exist in a 2×2 grid, so bending constraint list is empty here).
- Vertex 0 pinned. GUI radio buttons choose controlled vertex ∈ {1, 2, 3} (default 3).
- Loop (`:120-142`): `applyControls` **only while control is active** → `integrator->step(cloth, dt)` with **dt = 0.01, one step per frame** → **ad-hoc ground clamp:** any vertex with `y < 0.011` is set to `y = 0.01` and its y-velocity multiplied by 0.3 (this is *in addition to* the integrator's own y<0 clamp).
- GUI (ImGui): target X/Z sliders ∈ [−0.3, 0.5], Y ∈ [0.02, 0.8] (defaults 0.2, 0.3, 0.2); "Control Strength" slider ∈ [200, 2000] (default 800) used as the position gain; START calls `removeAllControls()` then `addPositionControl(controlVertex, target, strength, 100.0)` — **maxForce hardcoded to 100**; STOP; "Reset Cloth" (re-places at y=0.05, zeroes velocities, re-pins 0); "Drop Cloth" (places at y=0.3); live physics sliders: stiffness [100, 1500], damping [0.3, 0.95].
- **Master-branch difference (verified via `git diff master...tiny_opt -- optimization.cpp`):** on `master` this demo is `createGrid(10, 10, 0.1)` at y=0.5 with stiffness 800 / bending 20 / damping 0.9, controlling the whole `j=9` edge (10 vertices, indices `getGridIndex(i, 9)`) with directional pull buttons. **[README≠code]** The README's "Optimization" description (10×10 grid, edge pulling, pull-direction radio buttons) describes **master**, not this branch.

## 8. Build system

`CMakeLists.txt` — in-source builds are blocked; three executables, each compiled from its demo file **plus all of `src/*.cpp`** via `add_main_executable` (`CMakeLists.txt:148-194`): **`Visualization`**, **`Simulation`**, **`Optimization`** (output to `build/bin/`).

Actual dependency list:

| Dependency | How obtained | Required? | In README's install list? |
|---|---|---|---|
| Eigen3 ≥ 3.4 | system `find_package` (`:65`); vendored `deps/eigen` also on include path | REQUIRED | yes |
| TBB (oneTBB) | system, three-stage lookup (`:68-105`), CMake treats as *optional* | **effectively REQUIRED** — `src/integrator.cpp` calls `tbb::` unconditionally in `handleCollisions`; without `HAVE_TBB` the include is skipped and compilation fails | yes |
| OpenVDB | **system** `find_package(OpenVDB REQUIRED)` with hardcoded module path `/usr/local/lib/cmake/OpenVDB` (`:118-119`); vendored copy unused | REQUIRED (link-only; unused by code) | **no — README omits it** |
| LibXml2 | system `find_package(LibXml2 REQUIRED)` (`:145`) | REQUIRED (include-dir only; unused by code) | **no — README omits it** |
| Polyscope | `add_subdirectory(deps/polyscope)` | REQUIRED | (bundled) — needs OpenGL/GLFW/X11 dev packages, which README does list |
| Imath, Alembic | `add_subdirectory` (`:126`, `:129-135`), Alembic static, HDF5 off | built & linked, unused by code | (bundled) |
| OpenMP | optional fallback define only (`HAVE_OPENMP`); never actually used since the TBB code path is unconditional | optional | — |

Build steps (Linux/README): `mkdir build && cd build && cmake .. && make -j$(nproc)` → binaries in `build/bin/`. Other notes: `cmake_minimum_required(VERSION 3.0.0)` (`:1`) though README says 3.16+ **[README≠code]**; MSVC path exists in flags (`:41-55`) but the hardcoded `/usr/local/lib/cmake/OpenVDB` path and README target Ubuntu — **[inference]** Windows build likely needs manual OpenVDB/LibXml2 setup. `OpenCASCADE_INCLUDE_DIR` (`:163`) is referenced but never found — expands empty, harmless. No build directory exists in this checkout; nothing was verified by running.

## 9. Public API quick reference

### `SurfaceMesh` (base of `ClothMesh`)

| Method | Behavior |
|---|---|
| `void clear()` | drop all vertices/triangles |
| `void addVertex(const Vector3d&)` / `void addTriangle(int,int,int)` | append |
| `size_t getVertexCount() / getTriangleCount()` | counts |
| `const Vertex& getVertex(size_t)` | **no bounds check** |
| `void setVertexPosition(size_t, const Vector3d&)` | silently ignores bad index |
| `void computeNormals()` | area-weighted-ish vertex normals |
| `MatrixXd getVertexMatrix()` / `MatrixXi getTriangleMatrix()` | N×3 copies for Polyscope |

### `ClothMesh`

| Method / field | Behavior |
|---|---|
| `void createGrid(int width, int height, double spacing = 0.1)` | build grid at y=1.0, x=j·s, z=i·s; init state; build constraints |
| `int getGridIndex(int i, int j)` | `i * width + j` (row, col) |
| `void pinVertex(int)` | pin + zero velocity (no unpin API; use `pinned[i]=false`) |
| `void pinCorners()` | pin the 4 grid corners |
| `void addSphere(center, radius)` / `clearSpheres()` | collision spheres |
| `velocities, forces, pinned` | public per-vertex state vectors |
| `distanceConstraints, bendingConstraints` | public `vector<Edge>` |
| `properties` | `ClothProperties` (see §5 table) |

### `ClothController` — see §6 for signatures; one-liners

| Method | One-liner |
|---|---|
| `addPositionControl(v, target, gain=1000, maxForce=100)` | P-control toward target (creates/overwrites gains) |
| `addVelocityControl(v, vel, gain=100, maxForce=50)` | P-control on velocity; **overwrites shared maxForce** |
| `addForceControl(v, force)` | constant force each frame |
| `update{Position,Velocity,Force}Target(v, x)` | change setpoint only; no-op if not yet controlled |
| `applyControls(mesh, dt)` | advance controller clock; compute clamped F; `velocities[v] += F*dt` |
| `setTrajectory(v, positions, times, loop=false)` | linear-interp waypoint following; **resets gain/maxForce to 1000/100** |
| `addCircularMotion(v, center, r, f, axis=+Y)` / `addSinusoidalMotion(v, center, amp, f)` | targetPosition driven by pattern |
| `addWindForce(verts, dir, strength, turbulence=0)` | constant per-vertex force control |
| `removeControl(v)` / `removeAllControls()` | erase target + trajectory + pattern |
| `isControlled / getControlCount / getControlledVertices` | bookkeeping |
| `getControlForce(v)` | **only the externalForce component** (not PD force) |

### `SemiImplicitEulerIntegrator`

| Method | Behavior |
|---|---|
| `void step(ClothMesh&, double dt)` | full step: forces → integrate → collisions → 1 PBD iteration |
| `enableDebug(bool)` / `setDebugFrequency(int)` / `enableVerboseDebug(bool)` | console diagnostics |
| `const DebugInfo& getLastDebugInfo()` | totals from last debug computation |

## 10. Notes for the folding task

Facts most relevant to implementing symmetric and diagonal folds with this controller:

1. **Indexing/coordinates:** vertex `(row i → z, col j → x)` = index `i*W + j`; sheet lies in the xz-plane, y up, ground at y=0. For a W×H grid with spacing s the sheet spans `x∈[0,(W−1)s]`, `z∈[0,(H−1)s]`. Edge vertex sets: row `i=k` → indices `k*W + j (j=0..W−1)`; column `j=k` → `i*W + k (i=0..H−1)`.
2. **Set gravity explicitly** to `(0, −9.81, 0)` in any new demo — the struct default is (0,0,−9.81), which is horizontal. Also note effective gravity acceleration is `9.81/N` (§5), so free-fall phases of a fold are very slow on fine grids; plan trajectories that carry vertices down actively rather than relying on gravity to drop the folded flap.
3. **Driving vertices along a path:** for each grabbed vertex call `setTrajectory(v, waypoints, times, /*loop=*/false)` with `times[0] = 0.0`, then — because `setTrajectory` silently resets gains — call `addPositionControl(v, waypoints[0], yourGain, yourMaxForce)` **afterwards** to install your gains. Then per frame: `controller.applyControls(cloth, dt); integrator->step(cloth, dt);`. Trajectory time advances only via `applyControls` calls (controller-internal clock).
4. **Tracking/force budget:** pull force saturates at `maxForce` once position error exceeds `maxForce/gain` (defaults: 0.1 m). Per-frame velocity injection ≤ `maxForce*dt`. If waypoints move faster than the controller can track, the vertex lags and the effective path cuts corners **[inference — magnitude needs a run]**.
5. **No self-collision** (§5): the folded flap will pass through the stationary half without contact; the final "stacked" state is only shaped by the trajectory targets, the ground clamp at y=0 (integrator) and springs. Keep target heights ≥ one spacing above the resting half if you want visual layering; nothing enforces it physically.
6. **No shear constraints + only 1 PBD iteration on structural edges:** expect visible stretching and in-plane shearing during aggressive pulls, worst for the **diagonal fold** (diagonal motion is exactly the unresisted deformation mode). Countermeasures available in-code: raise `stiffness`, add substeps (Simulation demo uses 8×dt=0.003), or grab more vertices along the fold edge.
7. **Pinned vertices can't be moved by control** (velocity accumulates, position frozen — §6). For a fold, either pin nothing and control the anchor edge too, or unpin via `mesh.pinned[i] = false` before handing a vertex to the controller.
8. **Ground:** integrator clamps y<0 with no friction and no velocity zeroing; sliding is unresisted. The Optimization demo adds its own clamp at y=0.01 with y-velocity damping ×0.3 in the frame loop — replicate that pattern if you need calmer ground contact.
9. **This branch's Optimization demo is only 2×2** — for meaningful folds either build on master's 10×10 setup (stiffness 800, bending 20, damping 0.9, y=0.5, dt=0.01) or create a new grid.
10. **`FoldingOptimizer` (README §"TODO Guide") is purely aspirational.** None of `FoldingOptimizer`, `TrajectoryPoint`, `FoldDirection`, `Corner` exist anywhere in the codebase (grep-verified; the only "Corner" hit is `pinCorners`). All of these types and the class itself must be defined from scratch. The natural mapping: a `TrajectoryPoint` list per grabbed vertex → `ClothController::setTrajectory`, which already handles time-parameterized linear waypoint following.

## 11. Open questions and runtime unknowns

Only confirmable by building and running:

- Actual tracking fidelity/lag of position control at given (gain, maxForce, dt) — e.g. whether gain 800 / maxForce 100 tracks a 5 s fold trajectory closely at dt=0.01.
- How much the cloth stretches/shears under fold tension with only 1 constraint-projection iteration and no shear springs (statically expected to be significant; magnitude unknown).
- Stability envelope of semi-implicit Euler for stiffness 500–2000 at dt = 0.003–0.01 (demos suggest these are stable; margins unknown).
- Whether the unsynchronized parallel vertex writes in `handleCollisions` cause visible artifacts or nondeterminism in practice.
- Interaction when a controlled vertex is simultaneously inside a collision sphere or at the ground clamp (control Δv vs. position projection fighting each other).
- Real perceived gravity: the `g/N` scaling means large grids fall extremely slowly — how objectionable this looks for folding demos.
- Whether the build works at all on Windows/MSVC given the hardcoded Linux OpenVDB module path; and whether vendored-vs-system Eigen header mixing causes issues.
- Behavior of NaN-producing edge cases (circular motion with x-parallel axis; duplicate trajectory times).

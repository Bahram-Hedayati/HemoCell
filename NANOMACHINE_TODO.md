# Nanomachine molecular communication TODO

Repository analysis: 2026-09-22. Step 1 implemented in `examples/vessel`;
remaining transport/sensing work is listed below. See `examples/vessel/README.md`.

## Confirmed scope: first extension of pipeflow

Extend the copied project in `examples/vessel` with exactly two mobile nanomachines in
addition to its RBCs and platelets: one transmitter (TX) and one receiver (RX).
Both initially have the same size and shape as an RBC and move with the plasma
through HemoCell's existing immersed-boundary coupling. Their initial positions
are separated; their separation may evolve with flow and deformation.

TX releases glucose with the user-specified diffusion coefficient
`D_glucose = 9.2e-10 m^2/s`. Glucose is advected by the simulated plasma velocity
and diffuses. RX measures the glucose concentration on its moving surface.
The deliverable includes surface-resolved concentrations and an area-weighted
surface-average time series. Future machine shape/size changes must remain
possible through separate machine configuration.

Proposed modeling defaults still to specify: RBC-like mechanics as well as
geometry; continuum glucose concentration; passive sensing without glucose
removal; and release distributed over the TX exterior surface. Dose, timing,
initial positions/separation, background glucose, and membrane permeability
remain configuration/model decisions. Fixed probes and cell-free runs are
validation fixtures only, not substitutes for this milestone.

## What the repository already provides

| Area | Evidence and implications |
| --- | --- |
| Scalar transport | `config/constant_defaults.h` defines `CEPAC_DESCRIPTOR` as `AdvectionDiffusionD3Q19Descriptor`. `core/hemoCellFields.cpp::createCEPACfield()` creates an optional scalar lattice and registers `LatticeToPassiveAdvDiff3D` to couple fluid velocity to it. Reuse and validate this path before introducing a second solver. |
| Time integration | `core/hemoCell.cpp::iterate()` advances fluid, then CEPAC, then particles. Release and receiver sampling need explicit placement relative to those operations. |
| Existing release experiment | `examples/pipeflow/pipeflow.cpp` enables CEPAC, sets `tau_CEPAC = 1.0`, adds an RBC-mechanics `SYN` type, and every 10 steps resets a 3x3x3 box around each local SYN center to density 1.05. This is a prototype, not a calibrated molecule source. |
| Initialization discrepancy | The pipeflow comment says the scalar starts at zero, but `initializeAtEquilibrium` receives 1.0. Establish the descriptor's physical scalar convention before interpreting 1.0 or subtracting a background. |
| Configuration and scale | `examples/pipeflow/config.xml` enables CEPAC and sets `dx = 5e-7 m`, `dt = 1e-7 s`. Initially reuse the RBC mesh/resolution; reassess resolution when machine size/shape changes. |
| Existing output | `io/FluidHdf5IO.cpp` exports CEPAC through the fluid writer. `io/FluidHdf5IO.hh::outputDensity()` uses `computeDensity()` and, in SI mode, multiplies by `df/dx^2`; this is not a defined molecular concentration conversion. Pipeflow currently disables SI output globally. |
| Cell positions | `helper/cellInfo.{h,cpp}` provides centers, cell types, and `centerLocal`; useful for tracking TX/RX, but surface sensing also requires their mesh vertices/triangles and distributed ownership. |
| Restart gap | `core/hemoCellFields.cpp::save()` and `load()` serialize fluid and particle fields, but not CEPAC. Existing checkpoint calls do not preserve molecular transport history. |
| Dependency and input gaps | Root `setup.sh` pins and patches Palabos, but `palabos/` is absent in this checkout. Pipeflow references `SYN`, but its directory has no `SYN.xml` or `SYN.pos`; `core/hemoCellField.cpp` loads a material XML named after each type. |

No build or simulation was run during this analysis. Palabos constructor,
descriptor, coupling, and boundary behavior must be checked against the pinned
dependency before treating the existing transport configuration as validated.

## Example adaptation decision

- **Application: extend `examples/vessel`, copied from `examples/pipeflow`.** Retain its tube, plasma
  forcing, RBCs, and platelets. Replace the temporary generic SYN release
  experiment with distinct TX and RX roles and complete input files. Reuse
  CEPAC for glucose after validating its physical conversion and transport.
- **Validation reference: `examples/simple`.** Reuse its simple channel setup
  for isolated transport tests without STL or blood-cell complications.
- **API reference: `cases/CEPAC`.** Consult its scalar boundary setup, but its
  commented-out diffusion conversion and solidification coupling do not provide
  a validated glucose source or receiver model.

## Ordered implementation checklist

### 1. Add two RBC-shaped, mobile machines to vessel

- [x] Register TX and RX as distinct cell types/roles using `RbcHighOrderModel`
  and `RBC_FROM_SPHERE` initially, following the existing SYN registration.
  Supply `TX.xml`, `RX.xml`, `TX.pos`, and `RX.pos` with exactly one machine each.
  Initially copy RBC geometry/material settings; preserve existing RBC/PLT inputs.
- [x] Give machines stable cell IDs and explicit transmitter/receiver roles.
  Assert exactly one of each globally, including after restart and MPI migration.
  Do not identify machines by transient local particle indices.
- [ ] Configure initial centers/orientations and non-overlapping separation,
  checking clearance from other cells and vessel walls. Record initial and
  evolving separation, with periodic wrapping and unwrapped trajectories handled
  consistently. A fixed separation is not imposed on freely flowing machines.
- [x] Reuse existing force spreading, fluid-velocity interpolation, mesh movement,
  constitutive mechanics, and cell communication. Verify both machines deform
  and advect like RBCs and interact with RBCs/platelets/walls under the selected
  pipeflow interaction settings; assess whether repulsion must be enabled.
- [x] Keep geometry/material configuration separate from release/sensing roles.
  Avoid hard-coded radius, triangle count, or spherical sampling formulas so
  later shape/size changes can reuse the same communication implementation.
- [ ] Specify release dose/rate and timing, background glucose, initial machine
  placement, and observation duration. Use number concentration in molecules/m^3
  initially, with a documented scalar reference `C_ref` and optional molar output.
- [ ] Treat sensing as passive unless absorption is explicitly added later.
  The continuum scalar represents mean concentration, not individual arrival noise.

### 2. Validate and configure CEPAC transport

- [ ] Restore the pinned, patched Palabos dependency and establish a reproducible
  baseline build using the repository's CMake, MPI, and HDF5 setup.
- [ ] Inspect the pinned `AdvectionDiffusionBGKdynamics` constructor: determine
  whether it takes relaxation frequency `omega` or relaxation time `tau`.
  `createCEPACfield()` currently passes `param::tau_CEPAC` directly; the value
  1.0 hides a possible reciprocal error. Test a non-unit relaxation value.
- [ ] Set `D_glucose = 9.2e-10 m^2/s` in configuration and implement conversion in
  `mechanics/constantConversion.{h,cpp}` and XML configuration. Verify the
  descriptor constants before applying `D_lbm = D_phys*dt/dx^2` and
  `tau = 0.5 + D_lbm/c_s^2`, with `omega = 1/tau` if required by the API.
  Reject invalid settings; the current global `tau_CEPAC` default is zero.
  At the present pipeflow `dx` and `dt`, `D_lbm = 3.68e-4`; if the verified
  descriptor has `c_s^2 = 1/3`, this gives `tau = 0.501104` and
  `omega ~= 1.99559`. This near-limit relaxation requires targeted validation;
  do not replace the specified physical diffusion with a convenient value.
- [ ] Check scalar density/background conventions with uniform-field and
  zero-concentration tests. Resolve the pipeflow initialization discrepancy.
- [ ] Verify the velocity coupling executes before scalar transport with the
  intended fluid velocity, including external-force treatment. Check the
  expected equation `dc/dt + div(u*c) = div(D*grad(c)) + source`.
- [ ] Evaluate accuracy/stability at the desired diffusion and velocity scales;
  establish usable grid/time-step limits and refine or reconsider the scalar
  solver if the required regime is poorly resolved.
- [ ] Validate impermeable vessel walls for the scalar independently of fluid
  no-slip boundaries; do not assume the existing scalar bounce-back call is
  sufficient without a mass/leakage test.
- [ ] Make periodicity consistent across fluid, scalar, and machines.
  `setSystemPeriodicity()` currently updates fluid and particles, not CEPAC.
  For the periodic pipe, bound the observation interval and check domain-length
  sensitivity to recirculation; add clean inlet/outlet scalar conditions if
  a single-pass communication channel is required.

### 3. Define plasma/membrane transport and conservative TX release

- [ ] Decide and document glucose permeability for TX, RX, RBCs, and platelets.
  HemoCell's fluid IBM coupling does not itself impose molecular membrane
  boundaries: CEPAC occupies the Eulerian lattice, including cell interiors.
  Do not describe unrestricted scalar transport as plasma-only transport.
  If interiors are excluded, implement conservative moving-interface/no-flux
  treatment and account for changing accessible volumes; permeability requires
  its own interface law. Resolve this before finalizing deposition and sensing.
- [ ] Specify the release location: proposed first model is uniform release
  over the instantaneous TX exterior surface into plasma. Confirm whether a
  localized surface patch is needed. The existing center-box reset is not an
  exterior release model and may deposit inside TX.
- [ ] Implement an additive source compatible with scalar dynamics and local
  advection. For a pulse of `N` molecules require
  `sum(delta_c_i*accessible_volume_i) = N`; for a rate `q`, add `q*dt` each step.
  Preserve existing concentration and normalize surface quadrature/deposition
  weights so total dose is independent of deformation and mesh resolution.
- [ ] Update source support from TX's current mesh every release step. Handle
  walls, other cell interiors according to the permeability model, periodic
  wrapping, and empty support explicitly; never silently lose dose.
- [ ] Deposit exactly once into owned lattice nodes across blocks/ranks. A
  `centerLocal` condition does not cover a surface spanning multiple blocks.
  Avoid rank-dependent collective calls and synchronize scalar halos as needed.
- [ ] Define the ordering and time labels for release, fluid/scalar advancement,
  membrane motion, and RX sensing. Mesh and scalar samples must represent the
  same declared physical time. Include warmup, t=0 pulses, release intervals,
  and the existing particle/material update time-scale separations.

### 4. Sense glucose on the moving RX surface

- [ ] Interpolate glucose to the instantaneous receiver membrane at vertices
  or triangle quadrature points. Inspect `core/immersedBoundaryMethod.h` for
  reusable interpolation conventions, but use scalar-specific access and current
  mesh coordinates rather than stale fluid interpolation weights.
- [ ] Define the measured value as the plasma-side surface concentration.
  Choose a consistent one-sided interpolation/reconstruction if interiors are
  excluded or concentrations jump across membranes. Validate outward normals,
  wall proximity, and near-contact with other cells. If an exterior offset is
  numerically necessary, document it and demonstrate convergence toward the
  surface; do not substitute a center probe or arbitrary volume average.
- [ ] Output local surface concentrations and compute
  `C_RX(t) = integral_surface(c_plasma dA) / instantaneous_surface_area`
  using triangle quadrature or vertex-associated areas. Do not average vertices
  without area weighting on a nonuniform/deforming mesh. Optionally report
  surface minimum/maximum. A surface concentration integral is not a molecule
  count or absorbed dose.
- [ ] Count each surface element once globally, excluding ghost duplicates.
  Reduce area-weighted sums across MPI ranks and test receivers spanning blocks
  and periodic seams; reconstruct seam-crossing triangles consistently.
- [ ] Write a receiver CSV with iteration, physical time, TX/RX IDs and centers,
  separation, RX area, surface-average/minimum/maximum concentration, and
  cumulative released dose. Make sensing cadence independent of full-field I/O.
- [ ] Export surface concentrations associated with RX mesh IDs/coordinates and
  glucose field snapshots with physical units. Add concentration-specific HDF5
  conversion rather than disabling SI output globally for the fluid.
- [ ] Track total glucose, injected dose, boundary losses, and any explicit sinks
  or membrane exchange. Use the balance to detect numerical losses.

### 5. Integrate inputs, documentation, and restart support

- [ ] Update `examples/vessel/vessel.cpp`, its XML/position inputs, and README
  with the TX/RX scenario and build/run instructions. Use the `vessel` target registered in `examples/CMakeLists.txt`.
- [ ] Add validated XML settings for glucose diffusion, background/reference
  concentration, release schedule/dose, surface sensing, and sample intervals.
  Remove hard-coded temporary SYN behavior and provide all required input files.
- [ ] Document the chosen membrane model, surface sampling method, physical units,
  expected receiver output, and the limitations of a continuum glucose model.
- [ ] Persist glucose populations, machine identities, release counters, and
  sampling progress alongside fluid/particle checkpoints. Restore scalar coupling,
  boundaries, and surface ownership without duplicate releases/output records.
- [ ] Audit scalar/machine ownership through block restructuring and load
  balancing (`helper/loadBalancer.cpp`). Explicitly disable unsupported options
  until their transport state migration is implemented and tested.
- [ ] Choose a sufficient physical run time and sampling cadence for the initial
  separation and flow. The current 2,000 steps at `dt = 1e-7 s` cover only
  `2e-4 s`; determine the required duration from transport scales and pilot runs
  rather than assuming the existing configuration captures reception.

### 6. Validation and completion criteria

- [ ] Add focused tests under `tests/` and transport benchmarks under
  `tests/validation/`, using the existing GoogleTest/CMake arrangement.
- [ ] Verify zero/background preservation and exact injected dose, including
  repeated pulses, clipped source kernels, and non-unit relaxation values.
- [ ] Diffusion-only benchmark: compare a known initial distribution with the
  diffusion solution and variance growth; check mass and grid/time convergence.
- [ ] Uniform-flow benchmark: verify translation by `u*t` together with diffusion.
  Compare surface samples and the area-weighted receiver signal with the
  analytical field evaluated/integrated over the same test surface.
- [ ] Test no-flux wall mass conservation and periodic wrapping/recirculation;
  validate inlet/outlet mass accounting if open boundaries are added.
- [ ] Compare receiver traces and total dose on one and multiple MPI ranks,
  including moving TX/RX surfaces crossing partition boundaries and periodic seams.
- [ ] Compare uninterrupted and checkpoint/restart runs across a release event.
- [ ] Verify uniform concentration is recovered on every RX surface point and
  in its area-weighted mean, including deformed/nonuniform meshes. Test a known
  spatially varying field, mesh refinement, and plasma-side interpolation.
- [ ] Validate the chosen membrane treatment for stationary and moving surfaces,
  including no-flux conservation or prescribed permeability as appropriate.
- [ ] Use cell-free and fixed-surface cases only to isolate numerical errors.
  The acceptance run must retain RBCs and platelets plus exactly two freely
  mobile RBC-shaped TX/RX machines, releasing glucose at the specified diffusion
  coefficient and sensing on RX's actual moving surface.
- [ ] Check positivity, finite values, glucose mass balance, grid/time/mesh
  convergence, release support, particle-update cadence, and sensing cadence.
  Resolve accuracy/stability at the required low lattice diffusivity.
- [ ] Choose quantitative tolerances from refinement studies and document them.
  Completion means a reproducible run produces a physically calibrated
  surface-concentration-versus-time trace for mobile RX after a known TX dose,
  together with surface-resolved values and trajectories, with validated advection,
  diffusion, conservation, parallel execution, and restart behavior.

Deferred beyond the first validated version: receptor kinetics, chemical
reactions/decay, multiple species, active propulsion, customized machine
shape/size, and stochastic molecule-by-molecule reception. Mobile RBC-shaped
machines and sensing on their resolved mesh surfaces are part of the first version.

# Vessel: two mobile RBC-shaped nanomachines

This is the project directory for the nanomachine extension. Step 1 adds one
transmitter (`TX`) and one receiver (`RX`), using the same RBC mesh construction,
size, constitutive model and fluid coupling as `RBC`. They freely translate and
deform in the plasma. TX now releases glucose uniformly over its moving exterior surface at
100 molecules/s, with diffusion coefficient 9.2e-10 m²/s and advection by the
plasma velocity. RBC, platelet, TX and RX interiors are excluded from glucose
transport. RX surface sensing remains a later step.

The executable is built from `vessel.cpp`. The copied `pipeflow.cpp` is retained
as reference and is not compiled by this target. `examples/pipeflow` is unchanged.

## Initial configuration

The copied dense RBC/platelet position files are retained unchanged. They
contain 236 RBC and 17 platelet candidates; with the current tube geometry the
loader retains 36 complete RBCs and 6 platelets. Candidates outside the usable
geometry are rejected by HemoCell. Two additional machines are placed around
these loaded cells without replacing any RBC or platelet:

| Type | Count | Initial center (micrometres) |
| --- | --- | --- |
| TX | 1 | (5, 16, 18.5) |
| RX | 1 | (45.5, 16, 10) |

Both start with zero rotation. The center displacement from TX to RX is
(40.5, 0, -8.5) micrometres. The x direction is periodic (51.5 micrometres with
this voxelized configuration), so the shortest periodic center distance is
about 13.90 micrometres. Positions and separation evolve freely in the flow.

TX.xml and RX.xml copy RBC.xml's material and mesh settings, including the
3.91-micrometre radius parameter. Each `.pos` file starts with its cell count,
followed by `x y z rotationX rotationY rotationZ` rows (micrometres and degrees).
The radius in material XML is in metres. The machine positions were selected
using conservative separating-plane checks against the loaded initial RBC/PLT
meshes and periodic images. Changing geometry or blood-cell input requires
rechecking clearance; the loader does not prevent cell-cell overlap.

The copied dense files are also preserved in `initial_states/original_pipeflow/`.
The executable checks that exactly one complete TX and RX survives loading, at
each field output, and at the end of the run. Existing HemoCell cell IDs are
included in CSV and mesh output.

`config.xml` controls grid spacing, fluid flow, time step and output cadence.
`enableCEPACfield` must remain 0 for this step. Cell-cell repulsion is enabled
using the existing `kRep`/`RepCutoff` settings. Both machines use RBC's existing
material and particle update intervals. Separate TX/RX material files permit
later size/material changes; changing the shape construction will also require
changing the registration in `vessel.cpp`.

## Build

Replace `/path/to/HemoCell` below with the path to your clone.
Run commands in a terminal. The repository requires a C++ compiler, CMake, MPI
and HDF5 development libraries. On Debian/Ubuntu these are typically provided by
`build-essential cmake libopenmpi-dev openmpi-bin libhdf5-dev`. Python 3 with
NumPy and h5py is used for visualization conversion.

```bash
cd /path/to/HemoCell
# Only if palabos/ is missing: download the pinned source and apply HemoCell patches.
bash setup.sh
cmake -S . -B build -DBUILD_TESTING=OFF
cmake --build build --target vessel --parallel 4
```

Skip `setup.sh` when `palabos/` already exists: it replaces that dependency
folder. Alternatively, after dependency setup, run `bash compile.sh` from this
example. This configures CMake even when a build directory already exists.

If CMake cannot find `mpi.h` on this Ubuntu installation although OpenMPI is
installed, configure with its explicit header directory:

```bash
cmake -S /path/to/HemoCell -B /path/to/HemoCell/build \
  -DBUILD_TESTING=OFF \
  -DMPI_C_HEADER_DIR=/usr/lib/x86_64-linux-gnu/openmpi/include \
  -DMPI_CXX_HEADER_DIR=/usr/lib/x86_64-linux-gnu/openmpi/include
```

## Run

The working directory matters: geometry, material and position files are read
relative to the current directory.

```bash
cd /path/to/HemoCell/examples/vessel
mpirun -np 1 ./vessel config.xml
```

For two MPI ranks, replace `-np 1` with `-np 2`. Start a fresh run using this
version's inputs; old SYN/pipeflow checkpoints do not contain these four types
in the same order.

The default 2,000 iterations at `dt = 1e-7 s` simulate 0.2 milliseconds. Field
and CSV snapshots are written every 100 iterations, plus the initial state.
Check the log for `TX=1 RX=1 RBC=36 PLT=6`. Output is placed in `output/` or a
numbered variant if it already exists; use the directory printed by the run.
A short duration can give modest visible displacement. Increase `tmax` in
`config.xml` to observe more travel, checking counts and stability as you do.

`output/csv/TX.<iteration>.csv` and `RX.<iteration>.csv` contain centers,
velocities, area, volume and cell IDs. Positions are metres, velocity m/s,
area m^2, volume m^3; physical time is the filename iteration multiplied by dt.
Coordinates wrap at the periodic boundary, so an apparent position jump is not
loss of the machine.

```bash
python3 check_motion.py output
```

This checks finite CSV values, a constant blood-cell population, exactly one TX/RX, and persistent machine
IDs at every saved CSV time, and reports first-to-last displacement. Use the
actual output folder name. It is intended for the supplied short run
before a machine completes a full periodic circuit.

## Longer transport run (10 milliseconds)

`config_long.xml` preserves the physical and numerical parameters of `config.xml`
and changes only duration and output settings:

| Setting | Value | Physical interval |
| --- | --- | --- |
| `tmax` | 100000 | 0.01 s total |
| `tmeas`, `tcsv` | 1000 | 0.1 ms between saved frames |
| `tcheckpoint` | 5000 | 0.5 ms between checkpoints |
| Output directory | `output_long` | numbered automatically if already present |

Run from the example directory (no rebuild is required for XML changes):

```bash
cd /home/jorge/Bahram/HemoCell/examples/vessel
mpirun -n 1 ./vessel config_long.xml
```

There are 101 field frames including t=0. Based on the short validation's
approximately 0.11 s/step, budget roughly 3 hours, potentially longer as geometry
changes. Unconverted outputs are estimated at about 0.8–1 GB. These are estimates,
not a completed long-run benchmark. The 10 ms configuration has been checked for
consistent settings but has not yet been validated through all 100,000 steps.

At 100 molecules/s the final expected dose is 1 molecule. The free-space
one-axis diffusion standard deviation `sqrt(2*D*t)` is about 4.3 micrometres at
10 ms; confinement and cells change the actual distribution. The nominal maximum
fluid speed is about 0.022 m/s, so traversal of the 51.5-micrometre periodic length
can take only a few milliseconds. Glucose recirculates and accumulates because
there is no outlet or sink. This run examines the existing periodic model, not
an isolated single-pass vessel. Longer duration does not replace spatial/time
refinement validation or introduce individual molecular arrivals.

After completion, using the actual output directory printed in the log:

```bash
python3 check_glucose.py output_long
cd output_long
python3 ../../../scripts/FluidHDF5.py
python3 ../../../scripts/CellHDF5toXMF.py RBC
python3 ../../../scripts/CellHDF5toXMF.py PLT
python3 ../../../scripts/CellHDF5toXMF.py TX
python3 ../../../scripts/CellHDF5toXMF.py RX
```

Open `output_long/glucose/Glucose.*.xmf` directly for glucose and the converted
cell series from `output_long/`. Keep the color range fixed during playback;
rescale it using the longer run's data rather than the short run's maximum.
Physical frame spacing is now 0.0001 s. Periodic wrapping can make machines
appear to jump across the domain; the short-run endpoint displacement checker
is not an unwrapped long-run trajectory analysis.

If interrupted after a checkpoint, resume from the example directory:

```bash
mpirun -n 1 ./vessel output_long/checkpoint/checkpoint.xml
```

Use the same MPI rank count and keep the complete checkpoint together. If the
solver reports an unresolved exterior gap or invalid membrane, preserve the log
and checkpoint for diagnosis; extending time alone does not resolve that geometry.

## Convert output and visualize in ParaView

```bash
cd /path/to/HemoCell/examples/vessel/output
python3 ../../../scripts/FluidHDF5.py
python3 ../../../scripts/CellHDF5toXMF.py RBC
python3 ../../../scripts/CellHDF5toXMF.py PLT
python3 ../../../scripts/CellHDF5toXMF.py TX
python3 ../../../scripts/CellHDF5toXMF.py RX
```

Replace `output` with the run's actual output directory. The generated `.xmf`
files reference data in `hdf5/`; keep the files together.

1. Open ParaView, choose **File > Open**, and select the `TX.*.xmf` time series.
   Select the legacy **XDMF Reader** if prompted, as recommended by this
   repository's visualization guide. Click **Apply**.
2. Repeat for `RX.*.xmf`, `RBC.*.xmf` and `PLT.*.xmf`. Show each as **Surface**.
   Choose **Solid Color**, for example green TX, blue RX, red RBCs and yellow
   platelets. Separate output series make their roles identifiable without
   relying on numeric cell-type IDs.
3. Reset the camera, then use **Play** or the time slider to follow motion.
   Use **Surface With Edges** on TX/RX to inspect their RBC-shaped meshes.
4. Optionally open `Fluid.*.xmf`, apply a **Slice** through the tube center,
   and color it by **Velocity** magnitude. Keep the fluid hidden or translucent
   when it obscures the cell surfaces.
5. Compare the first and last TX/RX CSV centers as a numerical check of motion;
   visualization alone may not reveal very small displacements.

### Glucose concentration

Glucose XMF files are written directly; no postprocessing is needed for them.
Open `output/glucose/Glucose.*.xmf` as a series, choose **XDMF Reader** (legacy),
and click **Apply**. Color by `Concentration_molecules_per_m3` or
`Concentration_mol_per_m3`. Use a **Slice** through TX, or a **Contour** at a
positive concentration, together with the TX/RX surface series. Move beyond the
initial frame (which has zero glucose), then rescale the color range to data.
`PlasmaMask` is 1 in accessible plasma and 0 in excluded interiors/walls; use
**Threshold** on this array if displaying only accessible regions. Visualization
interpolation near masked surfaces is not a membrane permeability measurement.

Cell and glucose file series use the same saved iteration order. Physical time
is `iteration * dt`, recorded in glucose HDF5 attributes and `mass_balance.csv`.
Molar concentration is in mol/m³ (divide by 1,000 for mol/L).

```bash
python3 check_glucose.py output
```

This verifies every saved concentration field is finite and nonnegative, has
zero concentration in excluded nodes, and integrates to `100 * time_s` molecules.
Use `--rate VALUE` if you change the release rate.

## Glucose model and configuration

```xml
<glucose>
    <diffusionCoefficient>9.2e-10</diffusionCoefficient> <!-- m²/s -->
    <releaseRate>100</releaseRate> <!-- total molecules/s -->
</glucose>
```

Initial glucose is zero. Continuous release starts at simulation time zero after
fluid warmup and lasts until `tmax`. Each step adds `releaseRate * dt` expected
molecules. The default 0.0002-second run releases **0.02 expected molecules**.
This is a deterministic mean-concentration model: it does not simulate individual
molecular arrivals or fractional physical molecules. At 100 molecules/s a mean
of one molecule is released per 0.01 seconds, much longer than this demonstration.

`glucoseGrid.cpp` implements conservative finite-volume advection–diffusion on
the fluid grid, with limited linear reconstruction, upwind advection, central
diffusion and two-stage time stepping. Automatic substeps enforce its positivity
bound. The dimensionless diffusion coefficient is `D * dt / dx² = 0.000368`.
The code uses the specified physical D without substituting a larger value.
General finite-volume background: [LeVeque's materials](https://www.clawpack.org/fvmhp_materials/).

Release is area-weighted over the instantaneous TX triangles. Each triangle's
share is deposited into nearby visible exterior plasma nodes, with normalized
weights. Thus deformation changes the local flux per unit area while the total
rate stays 100 molecules/s. Each triangle is counted once, including across MPI
partitions and the periodic seam. Sources cannot deposit into cell interiors.

The scalar grid is periodic along x and has no flux through vessel walls or any
RBC/PLT/TX/RX membrane. Closed meshes determine interior masks; links intersecting
a membrane are blocked even when both endpoints are exterior. After movement,
glucose in newly covered grid volumes is conservatively redistributed along the
old plasma connectivity, using membrane-checked exterior paths if narrow gaps
are disconnected on the six-neighbor grid. An unresolved gap causes an explicit
error rather than silent dose loss. Glucose is neither consumed nor absorbed.

At each iteration, release and scalar transport use the current mesh/velocity;
HemoCell then advances fluid and particles; the scalar geometry is updated and
newly covered amounts remapped. Output describes this updated geometry. The
existing particle/material cadences are retained. `mass_balance.csv` records
dose, plasma total, excluded amount, time, TX area and transport substeps.
`remapped_molecules` is cumulative redistribution and may count the same mass
more than once; it is not extra release.

### Accuracy and scope

This is a reference implementation for the vessel demonstration. Cell interiors
use whole grid-volume exclusion, not fractional cut-cell volumes. Surface release
is represented by one quadrature point per triangle and a grid-scale exterior
stencil. Moving-boundary remapping, voxelized narrow gaps and numerical advection
diffusion introduce discretization error. Conservation and zero interior values
do **not** establish accurate near-surface concentration: grid, time-step and
mesh refinement studies are required before quantitative communication results.
The RX surface signal is not yet implemented.

HemoCell runs on MPI ranks, while velocity and owned membrane vertices are gathered
to rank zero for glucose updates. This avoids duplicate sources but limits memory
and speed scaling. Load balancing and changed MPI decompositions on restart have
not been validated. Keep this example's static decomposition. The periodic pipe
recirculates glucose; it is not a single-pass vessel with glucose-free inflow.
Keep `enableCEPACfield=0`: the existing CEPAC path does not impose these moving
impermeable membrane boundaries.

### Checkpoints and numerical tests

Glucose is saved as `checkpoint/glucose.h5` alongside HemoCell's checkpoint.
To resume, edit `tmax` in the saved `checkpoint/checkpoint.xml` to the desired
larger final iteration and run, from `examples/vessel`:

```bash
mpirun -n 1 ./vessel output/checkpoint/checkpoint.xml
```

Use the actual output folder. Keep all checkpoint files together and preserve
the XML declaration. Resume with the same grid, diffusion, release rate and MPI
rank count. Old checkpoints without glucose state cannot resume this extension.
Glucose output resumes without repeating its release or CSV row at the restart
iteration. The checkpoint `.old` backups must be restored as a consistent set
if manually rolling back a failed checkpoint write.

The standalone numerical tests need no MPI launch:

```bash
cd /path/to/HemoCell
cmake --build build --target vessel_glucose_test --parallel 4
./build/tests/vessel_glucose_test
```

They check physical diffusion variance, advection–diffusion refinement, dose,
positivity, blocked membranes, conservative movement, thin barriers, exterior
remapping and a seam-crossing transmitter.

## Validation performed

- Glucose extension: completed all 2,000 steps with 0.02 expected molecules,
  zero excluded-node concentration and final mass-balance error below 5e-16
  molecules. All 21 fields pass `check_glucose.py`. The successful validation
  output in this checkout is `output_glucose_0/`.
- One/two-rank glucose concentrations at step 100 agree to floating-point
  precision. A step-50 checkpoint resumed to step 100 agrees with uninterrupted
  output, with no repeated glucose diagnostic rows.
- ParaView 6.1.1's legacy reader loaded all 21 glucose frames through `pvpython`.
  The standalone numerical tests pass, including geometric remapping.
- Built the `vessel` CMake target against this checkout's pinned, patched Palabos.
- Completed the supplied 2,000-step configuration on one MPI rank: 36 RBCs,
  6 platelets, one TX and one RX persisted through all 21 CSV snapshots.
  TX moved approximately 2.677 micrometres along x; RX moved 3.173 micrometres.
- Completed a 100-step two-rank smoke test; TX/RX centers agree with the serial
  run at step 100 to the precision stored in CSV.
- Initial separating-plane checks confirmed nonintersection of TX/RX with the
  loaded RBC/PLT meshes and each other, including adjacent periodic images.
- All saved HDF5 datasets were finite. Generated 21 XMF snapshots per type and
  for the fluid, and checked their HDF5 references. ParaView GUI rendering has
  not been tested in this session.

Use the output directory produced by your latest run. During validation, all
21 timesteps for TX, RX, RBC, PLT and Fluid were loaded successfully using
ParaView 6.1.1's legacy XDMF reader through `pvpython` after the metadata fix
below. Generated output and executables are not included in Git.

## ParaView closes on Apply / Wayland warning

Two separate issues were found on this installation:

- The cell converter previously described `(N, 1)` arrays such as `Cell Id` as
  `Matrix` attributes with incomplete dimensions. This reproduced
  `Floating-point exception` with ParaView 6.1.1 even using the legacy reader.
  `scripts/CellHDF5toXMF.py` now writes these as node-centered scalar attributes
  with dimensions `N 1`. A regression test covers scalar/vector/matrix shapes.
- The installed ParaView bundle has `libqxcb.so` but no Wayland platform plugin.
  Launch it with the X11/XWayland backend explicitly:

  ```bash
  QT_QPA_PLATFORM=xcb paraview
  ```

Qt documents this per-process override at
<https://doc.qt.io/qt-6/qpa.html>. It does not change the desktop session.

No simulation rerun is needed to repair existing output. Regenerate older
XMF files from their HDF5 data after updating the converter:

```bash
cd /path/to/HemoCell/examples/vessel
../../scripts/batchPostProcess.sh 1
```

The `1` makes this repository script remove existing generated XMF files before
recreating them; it leaves HDF5 simulation data intact. Without it, the converter
skips existing XMF files, preserving any old malformed metadata. Then launch
ParaView using `QT_QPA_PLATFORM=xcb`, open each type as a separate file series,
select **XDMF Reader** (legacy), and click **Apply**. Start with `TX.*.xmf`.

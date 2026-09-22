# Vessel: two mobile RBC-shaped nanomachines

This is the project directory for the nanomachine extension. Step 1 adds one
transmitter (`TX`) and one receiver (`RX`), using the same RBC mesh construction,
size, constitutive model and fluid coupling as `RBC`. They freely translate and
deform in the plasma. Their names designate future communication roles; glucose
release, diffusion and surface sensing are not implemented in this step.

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

No concentration visualization should be expected yet. Glucose transport and
receiver surface-concentration output belong to the next implementation step.

## Validation performed

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

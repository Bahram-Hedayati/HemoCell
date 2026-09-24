# Longer-run membrane diagnostics

## Result

Two distinct failure mechanisms were reproduced in the earlier version.
Corrections are now written, but **have not been built or run**, at the user's
request. The recorded diagnostic results below describe the pre-correction code.

## Corrections awaiting user validation

- `core/voxelWallContact.h` implements a swept nearest-node voxel contact
  constraint. At a solid voxel, it suppresses blocked displacement components,
  retains permitted tangential motion and keeps the vertex on the fluid side.
  The full segment is checked, including corner crossings, so a step cannot
  tunnel through a wall just because its endpoint is fluid.
- `HemoCellParticleField::advanceParticles()` uses that constraint when
  `HemoCellFields::preserveMembranesAtWalls` is enabled. Vessel enables it for all
  cell types. Other examples keep the legacy behavior. This required changes
  outside the vessel directory; include the core files when committing.
- Contact does not delete a vertex or change mesh connectivity. Invalid starting
  positions inside solid still fail explicitly. The log counts owned-vertex
  corrections since process startup; this diagnostic counter resets on restart.
- `surfaceRay.h` uses half-open projected-triangle coverage with point-relative
  long-double edge predicates. Shared edges are counted consistently without
  merging distinct nearby intersections. Odd-count errors remain active.
- Tests cover the actual failing vertex position, tangential motion, detachment,
  thin walls, corners, negative coordinates, an open periodic seam, shared
  triangle edges, close paired crossings, and the recorded RBC 179 mesh.

This wall treatment is a kinematic constraint against the voxelized wall, not
an elastic/lubrication wall-force model. It changes near-wall trajectories and
can affect deformation. It does not introduce a matched reaction force into the
fluid or guarantee that entire triangle faces remain outside a nonconvex voxel
boundary. Vertex preservation alone is not validation of near-wall mechanics.
Grid refinement and contact-frequency/deformation checks remain necessary.
Multi-rank contact/ghost consistency still requires runtime validation; begin
with the one-rank reproduction. Existing glucose impermeability and conservation
checks remain enabled.

### Please build and run these commands

```bash
cd /home/jorge/Bahram/HemoCell
cmake --build build --target vessel vessel_glucose_test --parallel 4
./build/tests/vessel_glucose_test
```

Proceed only if the build and tests succeed. Then run a fresh 6,500-step case
with full glucose (do not add the membrane-only option):

```bash
cd /home/jorge/Bahram/HemoCell/examples/vessel
mpirun -n 1 ./vessel config_diagnostic.xml
python3 check_glucose.py output_diagnostic
python3 check_motion.py output_diagnostic
```

Use the actual numbered output directory printed by HemoCell, since earlier
`output_diagnostic` folders may already exist. Expect all four cell populations
to persist and glucose dose 0.065 expected molecules at the final step. A
nonzero wall-contact correction count is expected; it should be inspected rather
than treated as a mass loss or a deleted cell. Repeat with tmeas/tcsv both set
to 100 in a copy of the configuration, then extend duration only after these
checks pass. These expectations have not been verified on the corrected code.

### 1. Original incomplete-membrane error: wall contact

With the original long-run output interval of 1,000 steps, the membrane-only
reproduction stopped at **iteration 5,451** (0.5451 ms):

```
Incomplete glucose membrane at iteration 5451: cell=72 type=RBC
owned_vertices=641/642; missing vertex IDs: 205
```

Vertex 205 is absent from *all* current local and ghost records, not merely
excluded by the glucose ownership filter. Its last gathered position at step
5,450 was, in lattice units:

```
x = 47.383146577908008
y = 37.499993903119496
z = 47.524176434652446
```

This is far from the periodic x seam. The fluid boundary field identifies
node (47,37,48) as fluid and the adjacent node (47,38,48) as solid. The vertex
was only about 6.10e-6 lattice units below the y=37.5 nearest-node transition.
`HemoCellParticleField::advanceParticles()` advances each vertex, rounds its
position to a lattice node, tags vertices on solid nodes, then removes them.
Together, the missing vertex, its last position and boundary map identify wall
contact/deletion as the cause, rather than periodic ownership loss. This is
numerical mesh damage, not a model of physical RBC membrane rupture.

Vessel enables cell-cell repulsion but does not enable HemoCell boundary-particle
repulsion. The initial minimum wall clearance applies at loading; it does not
prevent later membrane contact during deformation. A suitable wall-interaction
correction must be selected and validated; simply disabling the glucose check
would violate the requested impermeable membrane model.

Diagnostic evidence is in `output_diagnostic_membranes/glucose/`, with prefix
`membrane_failure.000000005451`. These generated files are ignored by Git.
The membrane-only run retains the HemoCell physics and output cadence but omits
passive glucose computation; it is a diagnosis, not a molecular-result run.

### 2. Ray-intersection tolerance error

A full glucose run with snapshots every 100 steps stopped shortly after step
6,000. Restarting its step-5,000 checkpoint reproduced a ray error at **6,018**:
RBC 179, x-directed ray with transverse grid indices y=19, z=36.
All 642 vertices were present. Recomputing intersections from the saved double-
precision coordinates produced two legitimate hits:

```
81.46641267349136   triangle 116
81.46641270909512   triangle 669
separation: 3.56037617e-8 lattice units
```

`glucoseGrid.cpp` merges intersections separated by less than 1e-7 lattice units.
Here a near-edge grazing ray has two distinct crossings closer than that
threshold; merging them produces an odd hit count and a false unclosed-membrane
error. The corrective work should distinguish true shared-edge duplicates from
valid paired crossings, with a grazing-ray regression test. Merely ignoring odd
counts is not a safe correction.

Evidence is in `output_diagnostic_replay100/glucose/`, with prefix
`membrane_failure.000000006018` (gathered mesh iteration 6,018).

### Output cadence matters in this framework

`HemoCell::writeOutput()` also applies repulsion, synchronizes particles, deletes
incomplete cells and recalculates constitutive forces. Consequently it is not a
pure read-only observer. A replay switching to 1,000-step output at step 5,000
completed step 6,500, while the 100-step replay failed at 6,018. This is not proof
that changing output cadence fixes the simulation; preserve cadence when
reproducing a case and assess this dependence before quantitative conclusions.

## Changes and validation in this diagnostic pass

- More specific errors: iteration, cell/type, counts and missing IDs, or cell
  and ray coordinates for intersection failure.
- Read-only failure snapshots: current owned/ghost vertices, last gathered
  complete meshes with their iteration, and triangle connectivity.
- `config_diagnostic.xml`: 6,500-step reproduction, original 1,000-step output.
- `--diagnose-membranes-only`: faster cell/membrane collection diagnosis;
  intentionally disables scalar evolution/output and checkpoint saving.
- Built the vessel target. Diagnostic instrumentation agrees with the earlier
  validated glucose output at step 100 (relative L1 difference about 3.9e-18).
- The full diagnostic run maintained dose and exclusion through step 6,000.
  The 1,000-step-cadence checkpoint replay passed exported glucose checks through
  step 6,500, with total expected dose 0.065 molecules.

## Reproduce

From the repository root, build with:

```bash
cmake --build build --target vessel --parallel 4
cd examples/vessel
mpirun -n 1 ./vessel config_diagnostic.xml --diagnose-membranes-only
```

Use the same configuration without the final option to run full glucose.
Output is written under `output_diagnostic` or a numbered variant. For the
ray-error case, copy this configuration and set both tmeas and tcsv to 100.
The existing short `config.xml` is unchanged. Keep diagnostic outputs and
checkpoints when investigating further; generated output is not part of Git.

## Remaining validation work

1. Build and validate the new wall constraint; inspect corrections and near-wall
   deformation while retaining closed RBC/PLT/TX/RX meshes.
2. Run the new intersection regressions and confirm mass/interior-exclusion
   checks for both output cadences.
3. Re-test beyond both failure points, then test the complete intended duration.
   Do not claim the longer model is validated from the 6,500-step replay alone.

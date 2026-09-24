# Grazing-ray regression fixture

`rbc179_step6018.mesh` contains RBC 179 from the diagnostic checkpoint replay
that failed at iteration 6018. Coordinates are lattice units. The first line
contains the vertex and triangle counts, followed by xyz vertices in vertex-ID
order and zero-based triangle connectivity. The input is simulation-generated
geometry from this repository's dense vessel example, not a third-party dataset.

The old 1e-7 intersection-merging tolerance collapsed two crossings about
3.56e-8 lattice units apart into one. The standalone test checks the offending
ray and performs full geometry classification on this closed mesh. The fixture
is kept with the source so tests do not depend on ignored simulation outputs.

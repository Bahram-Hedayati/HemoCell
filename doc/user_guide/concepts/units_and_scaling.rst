Units and scaling
=================

HemoCell works internally in **lattice (LBM) units**, but you specify your
problem in **physical (SI) units**. Understanding the mapping between the two is
the single most important thing for setting up a stable, meaningful simulation.

The two quantities that anchor the mapping are set in the ``<domain>`` block of
your ``config.xml`` (see :ref:`xml_files:Config.xml`):

* ``<dx>`` — the physical length of one lattice cell, in metres.
* ``<dt>`` — the physical duration of one LBM time step, in seconds.

Together with the fluid density ``<rhoP>`` and kinematic viscosity ``<nuP>``,
these fix the conversion factors between SI and lattice units for every other
quantity (velocity, force, viscosity, ...). HemoCell computes the derived
lattice-Boltzmann parameters for you through the ``Parameters`` helper (aliased
as ``param`` in case files), for example::

   param::lbm_pipe_parameters((*hemocell.cfg), 50);
   param::printParameters();

Always call ``param::printParameters()`` early in a case so you can check the
resulting lattice values before committing to a long run.

.. note::

   ``<dx>`` defaults to **0.5 µm** in the example cases. Positions in ``.pos``
   files are given in micrometres, so a cell centre listed as ``10 12 24`` lands
   at lattice node ``20 24 48`` when ``dx = 0.5 µm``. This frequently surprises
   new users — see the :ref:`FAQ <FAQ:Frequently Asked Questions>`.

Choosing ``dx`` and ``dt`` for stability
----------------------------------------

The lattice Boltzmann method is only stable and accurate within a limited range
of its internal parameters. Two rules of thumb apply to essentially every
HemoCell simulation:

* keep the **lattice velocity** low: ``u_lbm < 0.1`` (well below the lattice
  speed of sound), and
* keep the **relaxation time** above its singular point: ``tau > 0.5`` (values
  comfortably above ``0.5``, e.g. ``~1.0``, are safest).

Because ``tau`` is determined by ``nuP``, ``dx`` and ``dt``, you tune stability
primarily by adjusting the resolution ``dx`` and the time step ``dt``. A
negative ``<dt>`` tells HemoCell to set ``tau = 1`` and back-calculate the
corresponding time step, which is a convenient starting point.

If a run shows jiggling cells, travelling waves, exploding densities, or cells
being deleted, the LBM parameters are the first thing to check — see
:ref:`common_mistakes:Numerical instabilities`.

Time-scale separation
---------------------

The fluid time step required for LBM stability is often *much* smaller than the
time step the cell mechanics would need on their own. Recomputing the membrane
constitutive model and re-interpolating velocities every single fluid step would
therefore be wasteful. HemoCell exploits this with **time-scale separation**:
expensive steps are run only once every *N* fluid iterations.

The relevant controls, set from a case file, are:

.. list-table::
   :header-rows: 1
   :widths: 45 55

   * - Call
     - Effect
   * - ``HemoCell::setMaterialTimeScaleSeparation()``
     - Re-evaluate a cell type's mechanical model every *N* fluid steps
       (the most expensive step).
   * - ``HemoCell::setParticleVelocityUpdateTimeScaleSeparation()``
     - Interpolate the fluid velocity onto the particles every *N* fluid steps.
   * - ``HemoCell::setRepulsionTimeScaleSeperation()``
     - Update the cell–cell repulsion force every *N* fluid steps.
   * - ``HemoCell::setInteriorViscosityTimeScaleSeperation()``
     - Update the interior-viscosity correction every *N* steps (and the
       expensive full-grid ray-tracing pass less frequently still).

Typical values in the examples are around ``20`` for the material model and
``5`` for the velocity update, but the right separation depends on your ``dt``
and flow conditions.

.. admonition:: Rule of thumb
   :class: tip

   Larger separations run faster but introduce a longer delay before cells
   "feel" the fluid (and vice versa). If you increase a separation and the
   simulation becomes unstable or cells lag the flow, reduce it again.

Force limiting
--------------

To keep a simulation stable even when a cell is momentarily subjected to extreme
deformation, the maximum force the constitutive model may exert on the fluid is
capped at compile time by ``FORCE_LIMIT`` (default **50 pN per surface point**,
see :ref:`constants`). This trades a little physical fidelity for unconditional
stability of the fluid coupling; it does *not* make an otherwise wrong model
correct.

.. _concepts:

Concepts
========

This section explains the ideas behind HemoCell: how the fluid and the cells are
modelled and coupled, how physical quantities map onto the simulation's internal
(lattice) units, how a cell is represented and given its mechanical behaviour,
and how a simulation is distributed across many processors.

If you just want to run a simulation, start with :ref:`from_source` and the
:ref:`hemocell_cases:Example cases`. Come back here when you want to understand
*why* a parameter exists or *how* to choose its value.

.. toctree::
   :maxdepth: 1

   overview
   units_and_scaling
   cell_model
   parallelism

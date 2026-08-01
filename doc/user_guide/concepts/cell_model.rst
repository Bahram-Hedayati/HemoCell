The cell model
==============

In HemoCell every cell is a **closed triangulated surface** — a mesh of vertices
(material points) connected by edges into triangles. The mesh represents the cell
membrane (and, for red blood cells, the combined membrane + cytoskeleton). The
mechanical behaviour of the cell emerges from forces defined on this mesh, which
are recomputed from its instantaneous shape at every material update.

Building the mesh
-----------------

When you add a cell type with ``HemoCell::addCellType()``, the second
argument selects how the initial mesh is constructed:

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Construction method
     - Description
   * - ``RBC_FROM_SPHERE``
     - Generate a biconcave red-blood-cell shape from the analytical equation in
       the HemoCell paper, refining a sphere to at least ``minNumTriangles``
       triangles.
   * - ``ELLIPSOID_FROM_SPHERE``
     - Generate an ellipsoid, mainly used for the platelet (PLT) model.
   * - ``WBC_SPHERE``
     - Generate a sphere representing a (white-blood-cell-like) nucleated cell.
   * - ``MESH_FROM_STL``
     - Load the mesh vertices from an STL file referenced in the cell's
       ``CELL.xml``.
   * - ``STRING_FROM_VERTEXES``
     - Legacy, no longer used.

The construction method, the target resolution (``minNumTriangles``), and the
reference geometry (``radius``, ``Volume``) are read from the cell's
``CELL.xml`` (see :ref:`xml_files:CELL.xml and CELL.pos`).

The membrane forces
-------------------

The constitutive (mechanical) model evaluates a set of force contributions over
the mesh. Most HemoCell cell models share the same four core terms, scaled by
coefficients you set per cell type in ``CELL.xml``:

.. list-table::
   :header-rows: 1
   :widths: 18 22 60

   * - Coefficient
     - Force term
     - Resists / represents
   * - ``kLink``
     - Link force
     - Stretching of the membrane along mesh edges (the in-plane elastic
       response, derived from a worm-like-chain / persistence-length model).
   * - ``kArea``
     - Local area force
     - Change in the area of individual triangles (local area conservation).
   * - ``kBend``
     - Bending force
     - Bending between adjacent triangles, i.e. the membrane's resistance to
       being curved (set in units of ``kBT``).
   * - ``kVolume``
     - Volume force
     - Change of the enclosed cell volume (quasi-incompressibility).

The deformation modes are bounded by global constants
(``MaxCellVolumetricChange``, ``MaxCellSurfaceAreaChange``,
``MaxCellBendingAngle``, ``MaxCellPersistenceLength``) defined in
``constant_defaults.h`` (the :ref:`cellgroup <reference/mechanical_models:Mechanical (cell) models>`
constants). These were tuned for validated RBC behaviour and improved stability
under high shear; see :cite:`Zavodszky:2017` and :cite:`tarksalooyeh2019optimizing`.

The actual force law lives in each model's ``ParticleMechanics`` method (in the
``mechanics/`` directory). Different models add, remove, or modify terms — for
example the white-blood-cell model adds an inner rigid-core term, while the
platelet model uses a much higher maximum bending angle. See
:ref:`reference/mechanical_models:Mechanical (cell) models` for the catalogue,
and the :ref:`FAQ <FAQ:Frequently Asked Questions>` for how to implement your
own.

Optional membrane physics
-------------------------

Two optional effects can be enabled per cell type, at the cost of extra
computation (both are compile-time features, see :ref:`constants`):

* **Interior viscosity** — gives the fluid *inside* the cell a different
  viscosity from the surrounding plasma, set via ``enableInteriorViscosity`` and
  ``viscosityRatio`` in ``CELL.xml``. Requires the library to be built with
  ``INTERIOR_VISCOSITY``. See
  :ref:`cases/cellCollision_interior_viscosity:Colliding cells with interior viscosity`.
* **Solidification mechanics** — used to model thrombus (clot) formation by
  freezing cells in place under certain conditions, enabled with
  ``HemoCell::enableSolidifyMechanics()`` and the
  ``SOLIDIFY_MECHANICS`` build option.

.. note::

   ``eta_m`` (membrane viscosity) appears in the cell models and ``CELL.xml``
   but is currently not used by the standard models. It is retained for models
   that implement a viscous membrane response.

Cell–cell and cell–wall interaction
-----------------------------------

Because of the immersed-boundary coupling, two membranes already interact when
they are about half an IBM kernel-width apart, before the meshes visually touch.
On top of this, an explicit short-range **repulsion** force can be enabled
between cells, and between cells and solid boundaries, to prevent overlap at high
hematocrit. See :ref:`advanced_cases/repulsion:Repulsive forces`,
``HemoCell::setRepulsion()``, and
``HemoCell::enableBoundaryParticles()``.

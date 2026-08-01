Frequently Asked Questions
==========================

Q: How do I add a different constitutive model for RBCs?
--------------------------------------------------------

The current constitutive models for all the cells are implemented under the ``HemoCell/mechanics`` folder.
For each model look for the implementation ``.cpp`` and the corresponding ``.h`` header file.
For instance, for RBCs look into the ``rbcHighOrderModel`` files. The implementation contains the function 
``ParticleMechanics`` that is called across the surface of the cell, and describes the forces arising from the deformation of the discretized (triangulated) surface of the cell. To create a new model duplicate the implementation and the header files under a different name, and redefine the contents of ``ParticleMechanics`` to any other desired model. Any model should work as long as it can be formulated as a force response function of deformation.

Q: Can I include the internal structural mechanics of a cell in the constitutive model?
---------------------------------------------------------------------------------------

It is possible to define arbitrary forces on the surface of a capsule / cell. For an example on how to mimic the mechanics of a nucleated cell (with a stiff nucleus) lok at the mechanical model of ``wbcHighOrderModel``. It represents a generic single-nucleus white blood cell.


Q: My cells are not initialized at the right position
-----------------------------------------------------

This is possible due to the lattice spacing usually (unless defined otherwise)
being 0.5 µm while the positions in the ``.pos`` file are in µm. This means
that a cell center of ``10 12 24`` in a ``.pos`` file will have a center of
``20 24 48`` in the simulation.

Q: Can the cells lie out of the domain?
---------------------------------------

Yes, cells are bound inside the domain by their centerpoints, so parts of the
cell might reach over the domain. These cells are deleted during initialisation,
unless the system is run on a single core.


Q: What is the exact centerpoint of the domain, how big is my domain?
---------------------------------------------------------------------

When you specify a domain of for example 10 cells, the exact
middle will be *BETWEEN* cell 5 and 6. imagine that dx is 0.5µm then the domain
is 9.5µm long, and the middle is 4.75µm. periodic boundaries of course add 1 dx
back to the length.

Q: Compiling with the singularity images gives an error within Cmake
---------------------------------------------------------------------

We have had reports that the clang (or llvm) compiler can cause problems with
compiling the singularity image. Disabling it (removing it from the environment)
or using another system might work. We are looking into this issue.

Q: Compiling with singularity gives errors
------------------------------------------

Maybe you have an old ``hemocell/`` folder in the same directory and are trying to compile with a newer version. If this is the case you can rename ``hemocell/`` to ``hemocell_old`` and the newer singularity image should create a new folder 

Q: HemoCell seems to be leaking memory over time
------------------------------------------------

We have experienced similar issues when compiling and running with OpenMPI
version 2.0.x, older or newer versions generally seem to be ok. We have tried to
track down the issue with heaptrack and only found memory leaks within the MCA
part of OpenMPI.

Q: Including custom example code raises: ``undefined reference``
----------------------------------------------------------------

This error might be encountered when linking external ``*.h,*.cpp`` files in
addition to your example's main ``cpp``-file without referencing these
files in the example's ``CMakeLists.txt`` file. Please refer to section
:ref:`linking_external_code` for a brief description how to link these files.

Q: Are there any undocumented features?
---------------------------------------
At the moment, sadly several of the new features lack documentation. These are often *researchware*, tested only partially, often in a single scenario. If you want to implement any new feature, please reach out first to make sure there is no unnecessary duplication of effort.


Do you have a problem that is not listed here? Please look at the `Discussion <https://github.com/UvaCsl/HemoCell/discussions>`_ section at the github repository! If you don't find your answer there open a new thread!

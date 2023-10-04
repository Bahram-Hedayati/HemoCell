Profiling Hemocell performance
================================

If you are interested in the performance of hemocell you can enable the inbuilt performance profiler.
The profiler measures the execution time of a selection of components per process.

The profiler is enabled by default when initializing the Hemocell object.

Adding the following line of code prints the measured performance of the first process to stdout.

.. code-block:: c++

    hemo::global.statistics.printStatistics();

Use this line of code to write the full profiling measurements to file.

.. code-block:: c++

    hemo::global.statistics.outputStatistics();

The default filename is "filename.statistics" to change the name add/change the ``config.xml`` file.

.. code-block:: xml

    <parameters>
        <filename>filename</filename>
    </parameters>


When running with a large number of processes writing the profiling results to file becomes very slow.
To parallelize writing the file you can add the optional ``batch_size`` parameter to the ``outputsStatistics()`` call.

.. code-block:: c++

    hemo::global.statistics.outputStatistics(batch_size);

When setting the ``batch_size``, the profiler output will be writting to multiple smaller files, which speeds up the process.
The value of ``batch_size`` determins howmany processes write to each file.
Decreasing this number will create a larger number of smaller files, instead of a single large file.
When in doubt set ``batch_size`` to the total number of processes per node. This creates a single statistics file per node.

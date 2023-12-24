#!/bin/bash
if [[ "$0" == "$BASH_SOURCE" ]]; then
 echo "Script is a subshell, this wont work, source it instead!"
 exit 1
fi

module load cmake
module load cray-hdf5
module load cray-python
module list

export CC=cc
export CXX=CC

# Compile by first running: cmake ..
# 	from hemocell/build directory


# NOTES
#
# Comment out line 130 in CMakeLists.txt to avoid messy unused parameter warnings.
#
# Improve the performnce of HemoCell by adding the following options to srun:
#
# srun --distribution=block:block --hint=nomultithread "$example" config.xml
#
# These options ensure you get the correct pinning of processes to cores on a compute node, without these options the default process placement may lead to a drop in performance for your jobs on ARCHER1.

#!/bin/bash
#
host=$(hostname)
if [ $host == pyro ]
then
    export NETCDF_HOME=/opt/netcdf-4.3.0/gnu/
    export NUMPY_HOME=/home/ronan/my27python/lib/python2.7/site-packages/numpy/core/
    export PYTHON_HOME=/home/ronan/my27python/include/python2.7
    export FOREFIRE_HOME=/home/ronan/Src/ForeFire/
elif [ $host == ubu ]
then
    export NETCDF_HOME=/opt/netcdf-4.3.0/gnu/
    export NUMPY_HOME=/usr/lib/python2.7/dist-packages/numpy/core/
    export PYTHON_HOME=/usr/include/python2.7
    export FOREFIRE_HOME=/home/paugam/Src/ForeFire/
else
    echo "hostname is not referenced"
fi
scons

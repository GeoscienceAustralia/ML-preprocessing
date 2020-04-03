#!/bin/env bash

# Fei Zhang 2020-04-03

# To test the latest code, run this script in gadi:
# ./run_tests_gadi.sh 
# fxz547@gadi-login-01 /g/data/ha3/fxz547/Githubz/geo-wavelets (master) $ 

module purge
module load python3/3.7.4
module load gdal/3.0.2
module load openmpi/2.1.6-mt
module load hdf5/1.10.5p
export GDAL_DATA=/apps/gdal/3.0.2/share/gdal/
export PYTHONPATH=/apps/gdal/3.0.2/lib64:/apps/gdal/3.0.2/lib64/python3.7/site-packages
export PATH=$HOME/.local/bin:$PATH
export PYTHONPATH=/g/data/ha3/fxz547/Githubz/geo-wavelets:$PYTHONPATH
export LC_ALL=en_AU.UTF-8
export LANG=en_AU.UTF-8


python3 -m pytest



exit


#############################################################################################
# example output shown below


================================================================ test session starts =================================================================
platform linux -- Python 3.7.4, pytest-5.2.2, py-1.8.0, pluggy-0.13.0
rootdir: /g/data/ha3/fxz547/Githubz/geo-wavelets
collected 100 items                                                                                                                                  

tests/test_multiscale.py ....................................................................................................                  [100%]

================================================================== warnings summary ==================================================================
/apps/gdal/3.0.2/lib64/python3.7/site-packages/osgeo/__init__.py:8
  /apps/gdal/3.0.2/lib64/python3.7/site-packages/osgeo/__init__.py:8: DeprecationWarning: the imp module is deprecated in favour of importlib; see the module's documentation for alternative uses
    import imp

/apps/gdal/3.0.2/lib64/python3.7/site-packages/osgeo/gdal.py:107
  /apps/gdal/3.0.2/lib64/python3.7/site-packages/osgeo/gdal.py:107: DeprecationWarning: gdal.py was placed in a namespace, it is now available as osgeo.gdal
    DeprecationWarning)

-- Docs: https://docs.pytest.org/en/latest/warnings.html
==================================================== 100 passed, 2 warnings in 107.56s (0:01:47) ========================

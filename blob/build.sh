#!/bin/bash

###############################################################################
# libhilbert_curve.a
###############################################################################

ln -s ../hilbert_curve/hilbert.cpp
ln -s ../hilbert_curve/hilbert.hpp
ln -s ../hilbert_curve/hilbert_c.cpp
ln -s ../hilbert_curve/hilbert_c.hpp
ln -s ../hilbert_curve_python_module/hilbert_curve_py.cpp
ln -s ../numpy_boost/numpy_boost.hpp
ln -s ../numpy_boost/numpy_boost_python.hpp

ln -s ../pyMDFortran/src/md.f90

# a soft link to a soft link links to the original file :-)
ln -s ../pyMD/experiment.py
ln -s ../pyMD/fcc.py
ln -s ../pyMD/logger.py
ln -s ../pyMD/particlecontainer.py
ln -s ../pyMD/verletlist.py

###############################################################################
# soft links into synchronized eclipse neon projects
###############################################################################

# 1. clean
rm -rf  ./hilbert.o ./hilbert_c.o  ./hilbert.d ./hilbert_c.d  libhilbert_curve.a

# 2. compile original hilbert code	
g++ -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"hilbert_c.d" -MT"hilbert_c.o" -o "hilbert_c.o" "hilbert_c.cpp"
 
# 3. compile wrappers of hilbert code
g++ -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"hilbert.d" -MT"hilbert.o" -o "hilbert.o" "hilbert.cpp"

# 4. link 1 and 2 into library
ar -r  "libhilbert_curve.a"  ./hilbert.o ./hilbert_c.o   

###############################################################################
# hilbert_curve_python_module
###############################################################################

# 1. clean
rm -rf  ./hilbert_curve_py.o  ./hilbert_curve_py.d  pyHilbertCpp.so
 
# 2. compile hilbert_curve_py.cpp 
g++ -I/Users/etijskens/miniconda3/envs/PP2016/lib/python3.5/site-packages/numpy/core/include -I/Users/etijskens/miniconda3/include/python3.5m -I/Users/etijskens/miniconda3/envs/PP2016/include -I"/Users/etijskens/ews/hpc-tnt-1.2/hilbert_curve" -O3 -Wall -c -fmessage-length=0 -Wno-unused-local-typedef -MMD -MP -MF"hilbert_curve_py.d" -MT"hilbert_curve_py.o" -o "hilbert_curve_py.o" "hilbert_curve_py.cpp"
 
# 3. build python module as shared library 
g++ -L./ -L/Users/etijskens/miniconda3/envs/PP2016/lib -L/Users/etijskens/miniconda3/envs/PP2016/lib/python3.5 -dynamiclib -o "pyHilbertCpp.so"  hilbert.o hilbert_c.o  hilbert_curve_py.o -lboost_python -lpython3.5m -lhilbert_curve

###############################################################################
# pyMDFortran python module
###############################################################################

# 1. clean
rm -f bin/testMDFortran *.mod

# 2. build Fortran module md.mod
gfortran -O2 -g -c md.f90

# 3. f2py to make pyMDFortran python module
f2py -c md.f90 -m pyMDFortran


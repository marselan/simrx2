#!/bin/bash

# compile simrx
rm -f simrx2
cd ..
make clean
make
cp simrx2 tests
cd tests

# test Water relative error
rm -f image1
mpiexec -np $1 simrx2 parameters_h2o.xml
octave test_water.m

# test Silicon relative error
rm -f image1
mpiexec -np $1 simrx2 parameters_si.xml
octave test_silicon.m
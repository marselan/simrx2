#!/bin/bash

# compile simrx
rm -f simrx2
rm -f image1
rm -f output_water
rm -f output_silicon
cd ..
make clean
make
cp simrx2 tests
cd tests

# test Water relative error
rm -f image1
rm -f output_water
mpiexec -np $1 simrx2 parameters_h2o.xml > output_water
octave test_water.m

# test Silicon relative error
rm -f image1
mpiexec -np $1 simrx2 parameters_si.xml > output_silicon
octave test_silicon.m

rm -f image1
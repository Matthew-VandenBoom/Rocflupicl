#!/bin/bash

# this causes bash to actually exit the script when a command results in an error, instead of just continuing for no reason
set -e

# module purge
# module load python/3.12 gcc/14.2.0 openmpi/5.0.7
# module list
cd libpicl
make clean
make
cd .. 
make clean 

rm -f build_lib/*.f90
rm -f build_lib/*.d
rm -f build_lib/*.o

rm -f build_util/*/*.f90
rm -f build_util/*/*.d
rm -f build_util/*/*.o

make -B RFLU=1 PICL=1 SPEC=1 FOLDER=1 DEBUG=1 -j16
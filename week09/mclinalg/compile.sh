#!/usr/bin/env bash

echo "Compile Fortran code to Python libs"
cd Fortran
make 
make clean
cd ..

echo "Run the python script"
cd Python
mpiexec -n 4 python numexperimentsmpi.py


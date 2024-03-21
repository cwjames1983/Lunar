#!/bin/bash


#Creates new output directories that could not be pushed to github:
#Should create the directories once and not override anything...

mkdir -p Outputs;
mkdir -p Outputs/Phis;
mkdir -p Outputs/Thetas;
mkdir -p Outputs/Time;
mkdir -p Outputs/Frequencies;
mkdir -p Outputs/Statistics;

g++ -O2 -o rough.exe rough.cpp -lm "/pawsey/mwa_sles12sp4/devel/cascadelake/gcc/8.3.0/openmpi-ucx-gpu/4.0.3/fftw/3.3.8/lib/libfftw3.a"

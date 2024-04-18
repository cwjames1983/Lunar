#!/bin/bash -l

#Creates new output directories that could not be pushed to github:
#Should create the directories once and not override anything...

mkdir -p Outputs;
mkdir -p Outputs/Phis;
mkdir -p Outputs/Thetas;
mkdir -p Outputs/Time;
mkdir -p Outputs/Frequencies;
mkdir -p Outputs/Statistics;

g++ -O2 -o rough.exe rough.cpp -lm "/usr/local/lib/libfftw3.a"

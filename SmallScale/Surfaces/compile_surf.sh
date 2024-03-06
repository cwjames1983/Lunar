#!/bin/bash

#06-02-24 Pretty sure this only works for my device. You'll need to configure
#the first command to match where you've downloaded fftw3.
#-cmw


# compiles rough surface generator
g++ -g -o surf.exe surf.cpp -lm /usr/local/lib/libfftw3.a

# compiles trivial surface generator
g++ make_trivial_surf.cpp -o trivial_surf.exe

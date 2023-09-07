# This command line specifies how to compile the rough.cpp code
# The file "rough_header.cpp" also needs to point to the
# location of fftw3.h, and the "Utilities" directory.
g++ -O2 -o rough.exe rough.cpp -lm "/usr/local/lib/libfftw3.a"

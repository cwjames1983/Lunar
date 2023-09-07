
# compiles rough surface generator
g++ -g -o surf.exe surf.cpp -lm /usr/local/lib/libfftw3.a

# compiles trivial surface generator
g++ make_trivial_surf.cpp -o trivial_surf.exe

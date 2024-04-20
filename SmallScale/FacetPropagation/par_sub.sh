#!/bin/bash -login

#SBATCH --job-name=ModD_par
#SBATCH --account=mwavcs
#SBATCH --partition=workq
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cores-per-socket=4
#SBATCH --cpus-per-task=4
#SBATCH --time=00:05:00
#SBATCH --output=Outputs/Terminal/%x.out

srun g++ -fopenmp -O2 -o rough.exe rough.cpp -lm "/pawsey/mwa_sles12sp4/devel/cascadelake/gcc/8.3.0/openmpi-ucx-gpu/4.0.3/fftw/3.3.8/lib/libfftw3.a"


export OMP_NUM_THREADS=4

#Change the input file to what you need.
srun ./rough.exe < Inputs/ModD/modD.in

#!/bin/bash -login

#SBATCH --job-name=ModD_par
#SBATCH --account=mwavcs
#SBATCH --partition=workq
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cores-per-socket=10
#SBATCH --cpus-per-task=10
#SBATCH --time=00:30:00
#SBATCH --output=Outputs/Terminal/%x.out

srun g++ -fopenmp -O2 -o rough.exe rough.cpp -lm "/pawsey/mwa_sles12sp4/devel/cascadelake/gcc/8.3.0/openmpi-ucx-gpu/4.0.3/fftw/3.3.8/lib/libfftw3.a"


export OMP_NUM_THREADS=10

#Change the input file to what you need.
srun -n 1 -c ${OMP_NUM_THREADS} ./rough.exe < Inputs/ModD/modD_fdiv.in

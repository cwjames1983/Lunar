#!/bin/bash -login

#SBATCH --job-name=ModD
#SBATCH --account=mwavcs
#SBATCH --partition=workq
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:10:00
#SBATCH --output=Outputs/Terminal/%x.out


#Change the input file to what you need.

srun ./rough.exe < Inputs/ModD/modD.in

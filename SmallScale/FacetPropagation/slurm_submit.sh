#!/bin/bash -l
#
#SBATCH --job-name=roughtest
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:00:30
#SBATCH --output=Outputs/Gary/%x.out

srun ./rough.exe < Inputs/quick_64.in

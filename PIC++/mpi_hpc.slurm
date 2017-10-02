#!/bin/bash
#SBATCH --nodes=2
#SBATCH --tasks-per-node=28
#SBATCH --cpus-per-task=1
#SBATCH -p tornado
#SBATCH -J PIC3d4
#SBATCH -o pic-%j.out
#SBATCH -e pic-%j.err
#SBATCH -t 240:00:00
mpirun -np 56 /home/ipntsr/romansky/PIC3d4/a.out

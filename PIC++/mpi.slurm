#!/bin/bash
#SBATCH -p all
#SBATCH -N 2
#SBATCH -n 24
#SBATCH -J PIC3d5
#SBATCH -o pic-%j.out
#SBATCH -e pic-%j.err
#SBATCH -t 480:00:00
scontrol show hostnames
hostname
mpirun -np 24 /home/uikam/romansky/PIC3d5/a.out

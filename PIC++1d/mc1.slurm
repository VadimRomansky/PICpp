#!/bin/bash
#SBATCH -p x5660
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -J mc
#SBATCH -o mc-%j.out
#SBATCH -e mc-%j.err
#SBATCH -t 480:00:00
scontrol show hostnames
hostname
./a.out

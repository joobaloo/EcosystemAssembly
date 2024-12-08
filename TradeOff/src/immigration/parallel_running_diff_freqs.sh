#!/bin/bash
#SBATCH--time=0-00:30:00
#SBATCH--ntasks=20
#SBATCH--cpus-per-task=1
#SBATCH--mem=300G
#SBATCH--partition=large_336 

srun --ntasks=1 echo "I'm task 1"
srun --ntasks=1 echo "I'm task 2"
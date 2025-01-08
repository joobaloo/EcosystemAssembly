#!/bin/bash
#SBATCH--time=0-05:00:00
#SBATCH--ntasks=21
#SBATCH--cpus-per-task=1
#SBATCH--mem=300G
#SBATCH--partition=large_336 
#SBATCH --output=job_%A_%a.log

# Array of immigration rates
rates=(10 20 40 80 160 320 640)

# Run simulations in parallel
for i in {0..6}; do
    srun --exclusive -n 1 julia ./src/immigration/assemble.jl 20 1 ${rates[i]} 1 1 5 &
done
wait

# Run analysis tasks in parallel
for i in {0..6}; do
    srun --exclusive -n 1 julia ./src/immigration/analysis.jl 20 1 ${rates[i]} 1 1 5 &
done
wait

# Run averaging in parallel
for i in {0..6}; do
    srun --exclusive -n 1 julia ./src/immigration/averages.jl 20 1 ${rates[i]} 1 1 5 &
done
wait

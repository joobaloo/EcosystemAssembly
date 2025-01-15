#!/bin/bash
#SBATCH --time=0-05:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=50G
#SBATCH --partition=large_336
#SBATCH --output=job_%A_%a.log
#SBATCH --array=0-6

# Array of immigration rates
rates=(10 20 40 80 160 320 640)

# Retrieve the rate based on the array task ID
rate=${rates[$SLURM_ARRAY_TASK_ID]}

# Determine the task type based on a passed argument
task_type=$1

if [[ "$task_type" == "simulation" ]]; then
    julia ./src/immigration/assemble.jl 20 1 $rate 1 20 25
elif [[ "$task_type" == "analysis" ]]; then
    julia ./src/immigration/analysis.jl 20 1 $rate 1 20 25
elif [[ "$task_type" == "averaging" ]]; then
    julia ./src/immigration/averages.jl 20 1 $rate 1 20 25
else
    echo "Invalid task type: $task_type"
    exit 1
fi


# Use these commands to run
#sbatch --array=0-6 src/immigration/array_parallel.sh simulation
#sbatch --array=0-6 src/immigration/array_parallel.sh analysis
#sbatch --array=0-6 src/immigration/array_parallel.sh averaging


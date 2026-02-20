#!/bin/bash
#SBATCH --time=0-10:00:00
#SBATCH --mem=70G
#SBATCH --partition=compute
#SBATCH --output=job_%A_%a.log
#SBATCH --array=0-2

cd /lustrehome/home/s.2540740/EcosystemAssembly/TradeOff

# Print debug info
echo "Running on: $(hostname)"
echo "Working directory: $(pwd)"


# Constant values
arg1=20
arg2=1
arg4=1
arg5=1
arg6=5

# # Values for arg3
# arg3_values=(10 80 320)

# # Select arg3 using SLURM_ARRAY_TASK_ID
# arg3=${arg3_values[$SLURM_ARRAY_TASK_ID]}
arg3=10
echo "Running with arguments: $arg1 $arg2 $arg3 $arg4 $arg5 $arg6"

julia --project=. src/immigration/analysis.jl "$arg1" "$arg2" "$arg3" "$arg4" "$arg5" "$arg6"
# max_retries=1
# retries=0

# while true; do
#     julia --project=. src/immigration/assemble_to_averages.jl \
#         "$arg1" "$arg2" "$arg3" "$arg4" "$arg5" "$arg6"

#     if [ $? -eq 0 ]; then
#         echo "Run succeeded."
#         break
#     else
#         retries=$((retries + 1))
#         echo "Error encountered. Attempt $retries of $max_retries."

#         if [ $retries -ge $max_retries ]; then
#             echo "Maximum retries reached."
#             break
#         fi

#         sleep 5
#     fi
# done
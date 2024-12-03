#!/bin/bash

# Set constant values
arg1=20       # Number of repeats
arg2=1        # Simulation type
arg4=1        # Number of immigrants
arg5=1        # Lower bound for number of reactions
arg6=25       # Upper bound for number of reactions

# Define the maximum number of retries
max_retries=3

# Define the array of values for the third argument (arg3)
arg3_values=(10 20 40 80 160 320)

# Loop over each value in the arg3_values array
for arg3 in "${arg3_values[@]}"; do
    retries=0  # Initialize retry counter

    echo "Running with arguments: $arg1 $arg2 $arg3 $arg4 $arg5 $arg6"

    # Retry loop
    while true; do
        # Run the Julia script with the current arguments
        julia ./src/immigration/assemble_to_averages.jl "$arg1" "$arg2" "$arg3" "$arg4" "$arg5" "$arg6"

        # Check the exit status of the Julia script
        if [ $? -eq 0 ]; then
            echo "Run succeeded with arguments: $arg1 $arg2 $arg3 $arg4 $arg5 $arg6"
            break
        else
            retries=$((retries + 1))
            echo "Error encountered. Attempt $retries of $max_retries."

            # Check if retries have reached the maximum allowed
            if [ $retries -ge $max_retries ]; then
                echo "Maximum retries reached. Moving to next set of arguments."
                break
            fi

            # Optional: Wait a few seconds before retrying
            sleep 5
        fi
    done
done

# End of script

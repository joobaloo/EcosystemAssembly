using JLD2
using Plots

function find_max_EUE_times()
    # Ensure at least one argument is passed
    if length(ARGS) < 1
        error("Need to provide an integer as the number of immigrants.")
    end

    # Parse the number of immigrants from the command line arguments
    num_immigrants = try
        parse(Int64, ARGS[1])
    catch e
        error("Need to provide a valid integer as the number of immigrants.")
    end

    println("Compiled and input read in!")
    flush(stdout)

    frequencies = [10, 20, 40, 80, 160, 320, 640]
    rl_vector = [1, 10, 20, 1]
    ru_vector = [5, 15, 25, 25]

    # Initialize an array to hold all max_EUE_times
    all_max_EUE_times = []

    # Open the JLD file and load the time data
    data_dir = joinpath(pwd(), "Output", "niche_size$(rl_vector[1])_$(ru_vector[1])", "1immigrants", "$(frequencies[1])events")

    stats_file = joinpath(data_dir, "RunStats$(frequencies[1])events_1immigrants.jld")

    # Check it actually exists
    if !isfile(stats_file)
        error("missing stats file for $(frequencies[1]) events 1 immigrant simulations")
    end

    # Extract time data
    t_times = load(stats_file, "times")

    for j in 1:4
        rl = rl_vector[j]
        ru = ru_vector[j]

        max_EUE_times = []

        for freq in frequencies
            # Construct the directory and file paths
            data_dir = joinpath(
                pwd(), "Output", "niche_size$(rl)_$(ru)", "$(num_immigrants)immigrants", "$(freq)events")
        
            stats_file = joinpath(data_dir, "RunStats$(freq)events_$(num_immigrants)immigrants.jld")

            # Check if the file exists
            if !isfile(stats_file)
                error("Missing stats file for $(freq) events, $(num_immigrants) immigrants simulations")
            end
        
            # Load simulation data
            t_times = load(stats_file, "times")
            community_EUE = load(stats_file, "mean_community_EUE")

            # Find the maximum EUE and corresponding time
            max_EUE_value, max_EUE_index = findmax(filter(!isnan, community_EUE))

            # Collect the time corresponding to the max EUE
            push!(max_EUE_times, t_times[max_EUE_index])
        end

        # Append the max_EUE_times for this (rl, ru) combination to the main array
        push!(all_max_EUE_times, max_EUE_times)
    end
    
    all_max_EUE_times_matrix = hcat(all_max_EUE_times)

    max_EUE_times_plot = plot(
        labels = ["1-5" "10-15" "20-25" "1-25"],
        xlabel = "Rate of Immigration",
        ylabel = "Time",
        title = "Time taken to reach max EUE"
    )
    for i in 1:4
        plot!(
            max_EUE_times_plot,
            frequencies,
            all_max_EUE_times[i],
            label = "niche_size$(rl_vector[i])_$(ru_vector[i])",
        )
    end
    # Define output directory and if necessary make it
    outdir = joinpath(pwd(), "Output", "Immigration_plots")
    mkpath(outdir)

    savefig(max_EUE_times_plot, joinpath(outdir, "max_EUE_times_plot_all_niche_sizes.png"))
    return (nothing)
end


@time find_max_EUE_times()
using Plots
using JLD2
using Statistics

function find_plateau()
     # Preallocate the variables I want to extract from the input
    num_immigrants = 0
    rl = 0
    ru = 0

    # Check that all arguments can be converted to integers
    try
        num_immigrants = parse(Int64, ARGS[1])
        rl = parse(Int64, ARGS[2])
        ru = parse(Int64, ARGS[3])
    catch e
        error("Need to provide an integer")
    end

    println("Compiled and input read in!")
    flush(stdout)

    frequencies = [10, 20, 40, 80, 160, 320, 640]

    # Open the JLD file and load the time data
    data_dir = joinpath(pwd(), "Output", "niche_size$(rl)_$(ru)", "1immigrants", "$(frequencies[1])events")

    stats_file = joinpath(data_dir, "RunStats$(frequencies[1])events_1immigrants.jld")

    # Check it actually exists
    if !isfile(stats_file)
        error("missing stats file for $(frequencies[1]) events 1 immigrant simulations")
    end

    # Extract time data
    t_times = load(stats_file, "times")

    rate_of_change_array = []
    threshold_values = []
    plateau_values_vector = []
    community_EUE_array = []

    for i in frequencies
        # Open the JLD file and load the surviving species data
        data_dir = joinpath(
        pwd(), "Output", "niche_size$(rl)_$(ru)", "$(num_immigrants)immigrants", "$(i)events")
    
        stats_file = joinpath(data_dir, "RunStats$(i)events_$(num_immigrants)immigrants.jld")

        # Check it actually exists
        if ~isfile(stats_file)
            error("missing stats file for $(i)events_$(num_immigrants)immigrants simulations")
        end
    
        # load simulation data
        t_times = load(stats_file, "times")
        community_EUE = load(stats_file, "mean_community_EUE")
       
        # differentiate to calculate rate of change in EUE
        rate_of_change = diff(community_EUE) ./ diff(t_times)

        # find plateau
        mean_roc = mean(rate_of_change)
        std_roc = std(rate_of_change)

        threshold = std_roc  # One standard deviation
        plateau_indices = findall(x -> abs(x - mean_roc) < threshold, rate_of_change)
        plateau_value = community_EUE[plateau_indices[1]]
        # collect data
        push!(rate_of_change_array, rate_of_change)
        push!(threshold_values, threshold)
        push!(community_EUE_array, community_EUE)
        push!(plateau_values_vector, plateau_value)

    end

    # Define output directory and if necessary make it
    outdir = joinpath(pwd(), "Output", "niche_size$(rl)_$(ru)", "Immigration_plots")
    mkpath(outdir)

     # initialise plots
    custom_layout = @layout [a b c d; e f g]
    rate_of_change_plot = plot(
        #legend = false,
        layout = custom_layout,
        plot_title="Rate of change of EUE, niche_size$(rl)_$(ru)",
        size=(2000,1000)
        )

    EUE_plot = plot(
        #legend = false,
        layout = custom_layout,
        plot_title="Plot plateau, niche_size$(rl)_$(ru)",
        size=(2000,1000)
        )

    # Add each frequency series to the plot with labels
    for (i, freq) in enumerate(frequencies)
        plot!(
            rate_of_change_plot, 
            t_times[1:end-1],
            rate_of_change_array[i], 
            xlims = (0, 6.3e7 / 2), 
            label="$(freq) events",
            seriestype = :line,
            subplot = i
        )

        hline!(
            rate_of_change_plot, 
            [threshold_values[i], -threshold_values[i]], 
            color=:red, 
            linestyle=:solid, 
            label="Threshold", 
            subplot = i
        )

        plot!(
            EUE_plot, 
            t_times, 
            xlims = (0, 6.3e7 / 2), 
            community_EUE_array[i],
            subplot = i, 
            label="$(freq) events"
        )

        hline!(
            EUE_plot,
            [plateau_values_vector[i]],
            color=:red, 
            linestyle=:solid,  
            subplot = i
        )
    end

    savefig(rate_of_change_plot, joinpath(outdir, " rate_of_change_plot.png"))
    savefig(EUE_plot, joinpath(outdir, "EUE_plateau_plot_$(frequencies[1])to$(frequencies[end])_frequencies.png"))
    
    return (nothing)
end

@time find_plateau()

using JLD2
using HypothesisTests
using Plots

using JLD
using Plots

function check_normality()
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
        error("Need to provide an integer: ", e)
    end

    println("Compiled and input read in!")
    flush(stdout)

    frequencies = [10, 20, 40, 80, 160, 320, 640]

    # Initialize array to collect data
    community_EUE_array = []

    for freq in frequencies
        # Construct the data directory and stats file paths
        data_dir = joinpath(pwd(), "Output", "niche_size$(rl)_$(ru)", "$(num_immigrants)immigrants", "$(freq)events")
        stats_file = joinpath(data_dir, "RunStats$(freq)events_$(num_immigrants)immigrants.jld")
        
        # Check if the file exists
        if !isfile(stats_file)
            error("Missing stats file for $(freq)events_$(num_immigrants)immigrants simulations")
        end

        # Load simulation data
        community_EUE = load(stats_file, "mean_community_EUE")

        # Collect data
        push!(community_EUE_array, community_EUE)
    end

    # Define output directory and if necessary make it
    outdir = joinpath(pwd(), "Output", "niche_size$(rl)_$(ru)", "Immigration_plots")
    mkpath(outdir)

    # # Initialize the main plot with a specific layout for subplots
    custom_layout = @layout [a b c; d e f g]
    distribution_plot = plot(
        layout= custom_layout, 
        plot_title="Distribution of EUE, niche_size$(rl)_$(ru)",
        size=(2000,1000)
        )

    # Loop to add subplots
    for (i, freq) in enumerate(frequencies)
        histogram_data = community_EUE_array[i]
        
        # Add a subplot to the corresponding position in the layout
        histogram!(distribution_plot, 
        histogram_data, 
        subplot = i, 
        title = "$freq",
        legend = false,
        xlims = (0, 1),
        xlabel = "EUE",
        ylabel = "Frequency"
        )
    end
    
    # Save the plot
    savefig(distribution_plot, joinpath(outdir, "distribution_of_EUE.png"))
end



function perform_anova()
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

    # Initialize arrays to collect data
    community_EUE_array = []

    for i in frequencies
        # Construct the data directory and stats file paths
        data_dir = joinpath(pwd(), "Output", "niche_size$(rl)_$(ru)", "$(num_immigrants)immigrants", "$(i)events")
        stats_file = joinpath(data_dir, "RunStats$(i)events_$(num_immigrants)immigrants.jld")
        
        # Check it actually exists
        if !isfile(stats_file)
            error("missing stats file for $(i)events_$(num_immigrants)immigrants simulations")
        end

        # Load simulation data
        community_EUE = load(stats_file, "mean_community_EUE")

        # Collect data
        push!(community_EUE_array, community_EUE)
    end

    # Combine all data into a single vector and create a corresponding group vector
    combined_data = vcat(community_EUE_array...)
    group_labels = repeat(1:length(frequencies), inner=length(community_EUE_array[1]))

    # Debug prints to verify data
    # println("Combined data length: ", length(combined_data))
    # println("Group labels length: ", length(group_labels))
    # println("Group labels: ", group_labels)

    # Perform the ANOVA test
    anova_result = OneWayANOVATest(combined_data, group_labels)

    # Inspect the results
    println(anova_result)
end

@time check_normality()
#@time perform_anova()

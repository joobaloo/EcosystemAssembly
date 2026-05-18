using JLD2
using Plots
using KernelDensity

function distribution_plots()
    # Preallocate the variables I want to extract from the input
    num_immigrants = 0

    # Check that all arguments can be converted to integers
    try
        num_immigrants = parse(Int64, ARGS[1])
        # rl = parse(Int64, ARGS[2])
        # ru = parse(Int64, ARGS[3])
    catch e
        error("Need to provide an integer: ", e)
    end

    println("Compiled and input read in!")
    flush(stdout)

    # Input vectors
    frequencies = [10, 80, 320]
    rl_vector = [1, 20, 1]
    ru_vector = [5, 25, 25]

    # Initialize a dictionary to collect data
    community_EUE_dict = Dict()

    for i in 1:3
        rl = rl_vector[i]
        ru = ru_vector[i]

        # Create a key for the dictionary
        data_key = (rl, ru)
        community_EUE_dict[data_key] = []

        for freq in frequencies
            # Construct the data directory and stats file paths
            data_dir = joinpath(pwd(), "Output", "niche_size$(rl)_$(ru)", "$(num_immigrants)immigrants", "$(freq)_a_year_rate")
            stats_file = joinpath(data_dir, "RunStats$(freq)_a_year_rate_$(num_immigrants)immigrants.jld")

            # Check if the file exists
            if !isfile(stats_file)
                error("Missing stats file for $(freq)_a_year_rate_$(num_immigrants) simulations")
            end

            # Load simulation data
            community_EUE = load(stats_file, "mean_community_EUE")

            # Collect data in the dictionary
            push!(community_EUE_dict[data_key], community_EUE)
        end
    end

    # Define output directory and if necessary make it
    outdir = joinpath(pwd(), "Output", "Immigration_plots")
    mkpath(outdir)

    # # Initialize the main plot with a specific layout for subplots
    custom_layout = @layout [a b; c d ; e f]
    distribution_plot = plot(
        layout= custom_layout, 
        plot_title="Distribution of EUEs",
        size=(2000,1000),
        bins=20, 
        )
    
    kde_plot = plot(
        layout= custom_layout, 
        #plot_title="Kernel Density Estimation of EUEs",
        size=(800,900),
        margin = 5Plots.mm
        )

    all_kde_plot = plot(
        layout= custom_layout, 
        plot_title="Kernel Density Estimation of EUEs",
        size=(2000,1000)
        )
    

    # Loop to add subplots
    for (i, freq) in enumerate(frequencies)
        EUE_data1_5 = community_EUE_dict[(1, 5)]
        #EUE_data10_15 = community_EUE_dict[(10, 15)]
        #EUE_data20_25 = community_EUE_dict[(20, 25)]
        #EUE_data1_25 = community_EUE_dict[(1, 25)]

        hist_data1_5 = EUE_data1_5[i]
        #hist_data10_15 = EUE_data10_15[i]
        #hist_data20_25 = EUE_data20_25[i]
        #hist_data1_25 = EUE_data1_25[i]

        kde_data1_5 = kde(hist_data1_5)
        #kde_data10_15 = kde(hist_data10_15)
        #kde_data20_25 = kde(hist_data20_25)
        #kde_data1_25 = kde(hist_data1_25)

        #all_kde = kde([hist_data1_5; hist_data20_25; hist_data1_25])
        #println(string(freq, " distribution spread is:", maximum(all_kde.x) - minimum(all_kde.x), " from ",  minimum(all_kde.x), " to ", maximum(all_kde.x)))
        
        
        # Add a subplot to the corresponding position in the layout

#

        # plot!(
        #     kde_plot,
        #     kde_data20_25.x, 
        #     kde_data20_25.density,
        #     label="Generalists",
        #     xlims=(-0.1, 1),
        #     #ylims=(0, 28),
        #     subplot = i,
        #     fillrange = 0,
        #     fillalpha = 0.4,
        #     color = :blue,
        #     legend = false
        # )

        # plot!(
        #     kde_plot,
        #     kde_data1_25.x, 
        #     kde_data1_25.density,
        #     label="Mix",
        #     xlims=(-0.1, 1),
        #     #ylims=(0, 28),
        #     tickfontsize = 12,
        #     subplot = i,
        #     fillrange = 0,
        #     fillalpha = 0.4,
        #     color = :red,
        #     legend = false
        # )

        # plot!(
        #     all_kde_plot,
        #     all_kde.x, 
        #     all_kde.density,
        #     label="all 3",
        #     xlims=(0, 1),
        #     ylims=(0, 22),
        #     subplot = i,
        #     fillrange = 0,
        #     fillalpha = 0.4,
        #     color = :red
        # )
        

    end
    
    # Save the plot
    #savefig(distribution_plot, joinpath(outdir, "grouped_distribution_of_EUE.png"))
    savefig(kde_plot, joinpath(outdir, "kde_of_EUE.png"))
    #savefig(all_kde_plot, joinpath(outdir, "all_niches_kde_of_EUE.png"))
end

function single_distribution_plots()
    # Preallocate the variables I want to extract from the input
    num_immigrants = 0
    rl = 0
    ru =0

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

    # This maps a Tuple of Ints to a Vector of Arrays
    community_EUE_dict = Dict{Tuple{Int, Int}, Vector{Any}}()

    # 2. Pre-seed the key so 'push!' has a target
    community_EUE_dict[(1, 5)] = []

    # Input vectors
    frequencies = [10, 80, 320]

    for freq in frequencies
            # Construct the data directory and stats file paths
            data_dir = joinpath(pwd(), "Output", "niche_size$(rl)_$(ru)", "$(num_immigrants)immigrants", "$(freq)_a_year_rate")
            stats_file = joinpath(data_dir, "RunStats$(freq)_a_year_rate_$(num_immigrants)immigrants.jld")

            # Check if the file exists
            if !isfile(stats_file)
                error("Missing stats file for $(freq)_a_year_rate_$(num_immigrants) simulations")
            end

            # Load simulation data
            community_EUE = load(stats_file, "mean_community_EUE")

            # Collect data in the dictionary
            push!(community_EUE_dict[(1,5)], community_EUE)
    end

    # Define output directory and if necessary make it
    outdir = joinpath(pwd(), "Output", "Immigration_plots")
    mkpath(outdir)

    # Create a 3-row layout for [10, 80, 320] frequencies
    custom_layout = @layout [a; b; c]
    kde_plot = plot(layout=custom_layout, size=(800, 900))

    for (i, freq) in enumerate(frequencies)
        # Access the pre-loaded data for the (1, 5) niche
        hist_data1_5 = community_EUE_dict[(1, 5)][i]
        
        # Generate the KDE curve
        kde_res = kde(hist_data1_5)

        plot!(
            kde_plot,
            kde_res.x, 
            kde_res.density,
            subplot = i,
            title = "Immigration Frequency: $freq",
            fillrange = 0,
            fillalpha = 0.4,
            color = :blue,
            legend = false,
            label = "Generalists (1, 5)",
            xlims = (-0.1, 1.1),
            ylims = (0, 25)
        )
    end

    savefig(kde_plot, joinpath(outdir, "kde_of_EUE.png"))

end

function single_distribution_histograms()
    # Preallocate variables
    num_immigrants = 0
    rl = 0
    ru = 0

    # Parse command line arguments
    try
        num_immigrants = parse(Int64, ARGS[1])
        rl = parse(Int64, ARGS[2])
        ru = parse(Int64, ARGS[3])
    catch e
        error("Need to provide an integer: ", e)
    end

    println("Compiled and input read in!")
    flush(stdout)

    # Use the dynamic key (rl, ru) as requested
    community_EUE_dict = Dict{Tuple{Int, Int}, Vector{Any}}()
    community_EUE_dict[(rl, ru)] = []

    frequencies = [10, 80, 320]

    for freq in frequencies
        data_dir = joinpath(pwd(), "Output", "niche_size$(rl)_$(ru)", "$(num_immigrants)immigrants", "$(freq)_a_year_rate")
        stats_file = joinpath(data_dir, "RunStats$(freq)_a_year_rate_$(num_immigrants)immigrants.jld")

        if !isfile(stats_file)
            error("Missing stats file for $(freq)_a_year_rate_$(num_immigrants) simulations")
        end

        # Load data
        community_EUE = load(stats_file, "mean_community_EUE")
        push!(community_EUE_dict[(rl, ru)], community_EUE)
    end

    outdir = joinpath(pwd(), "Output", "Immigration_plots")
    mkpath(outdir)

    # Layout for the 3 frequencies
    custom_layout = @layout [a; b; c]
    hist_plot = plot(layout=custom_layout, size=(800, 900))

    for (i, freq) in enumerate(frequencies)
        # Access data using the dynamic key
        raw_data = community_EUE_dict[(rl, ru)][i]
        
        # Plot raw distribution as a histogram
        histogram!(
            hist_plot,
            raw_data,
            subplot = i,
            bins = 30,             # Adjust bin count as needed for your data density
            title = "Immigration Frequency: $freq (Niche: $rl-$ru)",
            color = :green,
            alpha = 0.7,
            legend = false,
            xlabel = "Mean Community EUE",
            ylabel = "Frequency",
            xlims = (-0.1, 1.1),
            ylims = (0, 1200)
        )
    end

    savefig(hist_plot, joinpath(outdir, "histogram_of_EUE_$(rl)_$(ru).png"))
    println("Histogram saved to $(outdir)")
end


#@time distribution_plots()
#@time single_distribution_plots()
@time single_distribution_histograms()


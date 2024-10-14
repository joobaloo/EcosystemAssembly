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
    frequencies = [10, 20, 40, 80, 160, 320]
    rl_vector = [1, 10, 20, 1]
    ru_vector = [5, 15, 25, 25]

    # Initialize a dictionary to collect data
    community_EUE_dict = Dict()

    for i in 1:4
        rl = rl_vector[i]
        ru = ru_vector[i]

        # Create a key for the dictionary
        data_key = (rl, ru)
        community_EUE_dict[data_key] = []

        for freq in frequencies
            # Construct the data directory and stats file paths
            data_dir = joinpath(pwd(), "Output", "niche_size$(rl)_$(ru)", "$(num_immigrants)immigrants", "$(freq)events")
            stats_file = joinpath(data_dir, "RunStats$(freq)events_$(num_immigrants)immigrants.jld")

            # Check if the file exists
            if !isfile(stats_file)
                error("Missing stats file for $(freq)events_$(num_immigrants) simulations")
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
        EUE_data20_25 = community_EUE_dict[(20, 25)]
        EUE_data1_25 = community_EUE_dict[(1, 25)]

        hist_data1_5 = EUE_data1_5[i]
        #hist_data10_15 = EUE_data10_15[i]
        hist_data20_25 = EUE_data20_25[i]
        hist_data1_25 = EUE_data1_25[i]

        kde_data1_5 = kde(hist_data1_5)
        #kde_data10_15 = kde(hist_data10_15)
        kde_data20_25 = kde(hist_data20_25)
        kde_data1_25 = kde(hist_data1_25)

        all_kde = kde([hist_data1_5; hist_data20_25; hist_data1_25])
        println(string(freq, " distribution spread is:", maximum(all_kde.x) - minimum(all_kde.x), " from ",  minimum(all_kde.x), " to ", maximum(all_kde.x)))
        
        
        # Add a subplot to the corresponding position in the layout

        histogram!(
            distribution_plot, 
            hist_data1_5, 
            bins=10, 
            fillalpha = 0.5,
            label="1-5",
            title="$freq", 
            xlims=(0, 1),
            ylims=(0, 1800),
            xlabel="EUE", 
            ylabel="Frequency",
            subplot=i,
            color = :skyblue1
        )

        histogram!(
            distribution_plot, 
            hist_data20_25, 
            bins=10, 
            fillalpha = 0.5,
            label="20-25",
            title="$freq", 
            xlims=(0, 1),
            ylims=(0, 1800),
            xlabel="EUE", 
            ylabel="Frequency",
            subplot=i,
            color = :blue
        )

        histogram!(
            distribution_plot, 
            hist_data1_25, 
            bins=10, 
            fillalpha = 0.5,
            label="1-25",
            #title="$freq", 
            xlims=(0, 1),
            ylims=(0, 1800),
            xlabel="EUE", 
            ylabel="Frequency",
            subplot=i,
            color = :red
        )

        plot!(
            kde_plot,
            kde_data1_5.x, 
            kde_data1_5.density,
            label="Specialists",
            subplot = i,
            xlims=(-0.1, 1),
            #ylims=(0, 28),
            title="$freq/year", 
            fillrange = 0,
            fillalpha = 0.4,
            color = :skyblue1,
            legend = false
        )

        plot!(
            kde_plot,
            kde_data20_25.x, 
            kde_data20_25.density,
            label="Generalists",
            xlims=(-0.1, 1),
            #ylims=(0, 28),
            subplot = i,
            fillrange = 0,
            fillalpha = 0.4,
            color = :blue,
            legend = false
        )

        plot!(
            kde_plot,
            kde_data1_25.x, 
            kde_data1_25.density,
            label="Mix",
            xlims=(-0.1, 1),
            #ylims=(0, 28),
            tickfontsize = 12,
            subplot = i,
            fillrange = 0,
            fillalpha = 0.4,
            color = :red,
            legend = false
        )

        plot!(
            all_kde_plot,
            all_kde.x, 
            all_kde.density,
            label="all 3",
            xlims=(0, 1),
            ylims=(0, 22),
            subplot = i,
            fillrange = 0,
            fillalpha = 0.4,
            color = :red
        )
        

    end
    
    # Save the plot
    savefig(distribution_plot, joinpath(outdir, "grouped_distribution_of_EUE.png"))
    savefig(kde_plot, joinpath(outdir, "kde_of_EUE.png"))
    savefig(all_kde_plot, joinpath(outdir, "all_niches_kde_of_EUE.png"))
end

@time distribution_plots()



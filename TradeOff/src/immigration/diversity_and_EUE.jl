using JLD2
using Plots
using HypothesisTests
using StatsBase

function no_species_vs_EUE_plot()
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

    community_EUE_array = []
    num_species_array = []

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
        num_species = load(stats_file, "mean_surviving_species")

        # collect data
        push!(community_EUE_array, community_EUE)
        push!(num_species_array, num_species)

    end

    # Define output directory and if necessary make it
    outdir = joinpath(pwd(), "Output", "niche_size$(rl)_$(ru)", "Immigration_plots")
    mkpath(outdir)

    colour_palette = cgrad(:blues, length(frequencies))

    num_species_vs_EUE_plot = plot(
        xlabel = "no. species",
        ylabel = "EUE",
        ylims = (0, 1),
        title = "species richness and EUE, niche_size$(rl)_$(ru)",
        legend = true
    )

    # Add each frequency series to the plot with labels
    for (i, freq) in enumerate(frequencies)   
        plot!(
            num_species_vs_EUE_plot, 
            num_species_array[i], 
            community_EUE_array[i], 
            label="$(freq)/year",  
            ms=2, 
            ma=0.5,
            msw = 0,
            seriestype = :scatter,
            color = colour_palette[i],
            #legend = false
            )
    end

    savefig(num_species_vs_EUE_plot, joinpath(outdir, "num_species_vs_EUE_plot_$(frequencies[1])to$(frequencies[end])_frequencies.png"))
    
    # spearmans rank correlation
    community_EUE_vector = vcat(community_EUE_array...)
    num_species_vector = vcat(num_species_array...)
    correlation_results = corspearman(community_EUE_vector, num_species_vector)
    println(correlation_results)
    return (nothing) 
end

function shannon_diversity_vs_EUE_plot()
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

    frequencies = [10, 20, 40, 80, 160, 320]

    # Open the JLD file and load the time data
    data_dir = joinpath(pwd(), "Output", "niche_size$(rl)_$(ru)", "1immigrants", "$(frequencies[1])events")

    stats_file = joinpath(data_dir, "RunStats$(frequencies[1])events_1immigrants.jld")
    
    # Check it actually exists
    if !isfile(stats_file)
        error("missing stats file for $(frequencies[1]) events 1 immigrant simulations")
    end

    # Extract time data
    t_times = load(stats_file, "times")

    community_EUE_array = []
    shannon_diversity_array = []

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
        shannon_diversity = load(stats_file, "mean_shannon_diversity")

        # collect data
        push!(community_EUE_array, community_EUE)
        push!(shannon_diversity_array, shannon_diversity)

    end

    # Define output directory and if necessary make it
    outdir = joinpath(pwd(), "Output", "niche_size$(rl)_$(ru)", "Immigration_plots")
    mkpath(outdir)

    colour_palette = cgrad(:blues, length(frequencies))

    shannon_diversity_vs_EUE_plot = plot(
        xlabel = "Shannon Diversity Index (H)",
        ylabel = "EUE",
        ylims = (0, 1),
        #title = "shannon diversity and EUE, niche_size$(rl)_$(ru)",
        legend = false
    )

    # Add each frequency series to the plot with labels
    for (i, freq) in enumerate(frequencies)   
        plot!(
            shannon_diversity_vs_EUE_plot, 
            shannon_diversity_array[i], 
            community_EUE_array[i], 
            label="$(freq)/year",
            ms=2, 
            ma=0.7,
            msw = 0,
            tickfontsize = 12,
            seriestype = :scatter,
            color = colour_palette[i]
            )
    end

    savefig(shannon_diversity_vs_EUE_plot, joinpath(outdir, "shannon_diversity_vs_EUE_plot_$(frequencies[1])to$(frequencies[end])_frequencies.png"))
    
    # spearmans rank correlation
    community_EUE_vector = vcat(community_EUE_array...)
    shannon_diversity_vector = vcat(shannon_diversity_array...)
    correlation_results = CorrelationTest(community_EUE_vector, shannon_diversity_vector)
    println(correlation_results)
    correlation_results_v2 = corspearman(community_EUE_vector, shannon_diversity_vector)
    println(correlation_results_v2)
    return (nothing) 
end



#@time no_species_vs_EUE_plot()
@time shannon_diversity_vs_EUE_plot()
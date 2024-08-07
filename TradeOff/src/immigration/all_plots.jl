using JLD2
using Plots
using StatsPlots

function all_plots()
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
        num_substrates_array = []
        total_biomass_array = []

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
            num_substrates = load(stats_file, "mean_no_substrates")
            total_biomass = load(stats_file, "mean_total_biomass_of_viable_species")

            # collect data
            push!(community_EUE_array, community_EUE)
            push!(num_species_array, num_species)
            push!(num_substrates_array, num_substrates)
            push!(total_biomass_array, total_biomass)

        end

        # Define output directory and if necessary make it
        outdir = joinpath(pwd(), "Output", "niche_size$(rl)_$(ru)", "Immigration_plots")
        mkpath(outdir)
        
        # initialise plots
        EUE_plot = plot(xlabel="Time", ylabel="EUE")
        num_species_plot = plot( xlabel="Time",ylabel="Number of species")
        num_substrates_plot = plot( xlabel="Time",ylabel="Number of substrates")
        biomass_plot = plot( xlabel="Time",ylabel="Total Biomass")

        # Add each frequency series to the plot with labels
        for (i, freq) in enumerate(frequencies)
            plot!(EUE_plot, t_times, xlims = (0, 6.3e7 / 2), community_EUE_array[i], label="$(freq) events")
            plot!(num_species_plot, t_times, xlims = (0, 6.3e7 / 2), num_species_array[i], label="$(freq) events")
            plot!(num_substrates_plot, t_times, xlims = (0, 6.3e7 / 2), num_substrates_array[i], label="$(freq) events")
            plot!(biomass_plot, t_times, xlims = (0, 6.3e7 / 2), total_biomass_array[i], label="$(freq) events")
        end

        p = plot(
            EUE_plot, 
            num_species_plot, 
            num_substrates_plot, 
            biomass_plot, 
            layout = (2, 2), 
            plot_title = "niche_size:$(rl)-$(ru), $(frequencies[1])to$(frequencies[end])_rates, 1 immigrant, 10 repeats",
            size = (1200, 800),
            margin = 10Plots.mm
            )

        savefig(p, joinpath(outdir, "all_4_plots_$(frequencies[1])to$(frequencies[end])_frequencies.png"))
        savefig(EUE_plot, joinpath(outdir, "EUE_$(frequencies[1])to$(frequencies[end])_frequencies.png"))
        savefig(num_species_plot, joinpath(outdir, "num_species_$(frequencies[1])to$(frequencies[end])_frequencies.png"))
        savefig(num_substrates_plot, joinpath(outdir, "num_substrates_$(frequencies[1])to$(frequencies[end])_frequencies.png"))
        savefig(biomass_plot, joinpath(outdir, "biomass_$(frequencies[1])to$(frequencies[end])_frequencies.png"))
        return (nothing)
end

function final_max_EUE_plots()
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
        error("missing stats file for $(frequencies[1])events 1 immigrant simulations")
    end

    # Extract time data
    t_times = load(stats_file, "times")

    community_EUE_array = []
    max_EUE_times = []

    final_EUE = []
    final_EUE_SE =[]
    max_EUE = []
    max_EUE_SE =[]

    final_surviving_species = []
    final_surviving_species_SE = []
    max_surviving_species = []
    max_surviving_species_SE = []

    final_total_biomass = []
    final_total_biomass_SE = []
    max_total_biomass = []
    max_total_biomass_SE = []
    
    final_no_substrates = []
    final_no_substrates_SE = []
    max_no_substrates = []
    max_no_substrates_SE = []


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
        community_EUE_SD = load(stats_file, "sd_community_EUE")
        surviving_species = load(stats_file, "mean_surviving_species")
        surviving_species_SD = load(stats_file, "sd_surviving_species")
        total_biomass = load(stats_file, "mean_total_biomass_of_viable_species")
        total_biomass_SD = load(stats_file, "sd_total_biomass_of_viable_species")
        no_substrates = load(stats_file, "mean_no_substrates")
        no_substrates_SD = load(stats_file, "sd_no_substrates")
    
        max_EUE_value, max_EUE_index = findmax(filter(!isnan, community_EUE))
        max_surviving_species_value, max_surviving_species_index = findmax(filter(!isnan, surviving_species))
        max_total_biomass_value, max_total_biomass_index = findmax(filter(!isnan, total_biomass))
        max_no_substrates_value, max_no_substrates_index = findmax(filter(!isnan, no_substrates))

        # collect data
        push!(community_EUE_array, community_EUE)
        push!(max_EUE_times, t_times[max_EUE_index])

        push!(final_EUE, community_EUE[end])
        push!(final_EUE_SE, community_EUE_SD[end]/3)
        push!(max_EUE, max_EUE_value)
        push!(max_EUE_SE, community_EUE_SD[max_EUE_index]/3)
        
        push!(final_surviving_species, surviving_species[end])
        push!(final_surviving_species_SE, surviving_species_SD[end]/3)
        push!(max_surviving_species, max_surviving_species_value)
        push!(max_surviving_species_SE, surviving_species_SD[max_surviving_species_index]/3)

        push!(final_total_biomass, total_biomass[end])
        push!(final_total_biomass_SE, total_biomass_SD[end]/3)
        push!(max_total_biomass, max_total_biomass_value)
        push!(max_total_biomass_SE, total_biomass_SD[max_total_biomass_index]/3)

        push!(final_no_substrates, no_substrates[end])
        push!(final_no_substrates_SE, no_substrates_SD[end]/3)
        push!(max_no_substrates, max_no_substrates_value)
        push!(max_no_substrates_SE, no_substrates_SD[max_no_substrates_index]/3)
    end

    # Define output directory and if necessary make it
    outdir = joinpath(pwd(), "Output", "niche_size$(rl)_$(ru)", "Immigration_plots")
    mkpath(outdir)

    EUE_matrix = hcat(final_EUE, max_EUE)        
    EUE_SD_matrix = hcat(final_EUE_SE, max_EUE_SE)

    no_species_matrix = hcat(final_surviving_species, max_surviving_species)        
    no_species_SD_matrix = hcat(final_surviving_species_SE, max_surviving_species_SE)

    total_biomass_matrix = hcat(final_total_biomass, max_total_biomass)        
    total_biomass_SD_matrix = hcat(final_total_biomass_SE, max_total_biomass_SE)

    no_substrates_matrix = hcat(final_no_substrates, max_no_substrates)        
    no_substrates_SD_matrix = hcat(final_no_substrates_SE, max_no_substrates_SE)

    final_max_EUEs_bar_plot = groupedbar(
        EUE_matrix,
        yerr = EUE_SD_matrix,
        title = "Final and max EUEs, niche_size$(rl)_$(ru)", 
        xlabel = "Rate", 
        xticks = (1:length(frequencies), frequencies),
        ylabel = "EUE",
        ylims = (0, 1),
        label = ["final" "max"],
        bar_position = :dodge
        )
    
    final_max_EUEs_plot = plot(
        EUE_matrix,
        ribbon = EUE_SD_matrix,
        label = ["final" "max"],
        xlabel="Rates",
        xticks = (1:length(frequencies), frequencies),
        ylabel ="EUE",
        ylims = (0, 1),
        #title = "final and max EUE, niche_size$(rl)_$(ru)" 
        )
    
    final_max_no_species_plot = plot(
        no_species_matrix,
        ribbon = no_species_SD_matrix,
        label = ["final" "max"],
        xlabel="Rates",
        xticks = (1:length(frequencies), frequencies),
        ylabel ="no. species",
        #title = "final and max no. species, niche_size$(rl)_$(ru)" 
        )

    final_max_total_biomass_plot = plot(
        total_biomass_matrix,
        ribbon = total_biomass_SD_matrix,
        label = ["final" "max"],
        xlabel="Rates",
        xticks = (1:length(frequencies), frequencies),
        ylabel ="total biomass",
        #title = "final and max total biomass, niche_size$(rl)_$(ru)" 
        )
    
    final_max_no_substrates_plot = plot(
        no_substrates_matrix,
        ribbon = no_substrates_SD_matrix,
        label = ["final" "max"],
        xlabel="Rates",
        xticks = (1:length(frequencies), frequencies),
        ylabel ="no. substrates",
        #title = "final and max no. substrates, niche_size$(rl)_$(ru)" 
        )

    all_4_final_max_plots = plot(
        final_max_EUEs_plot,
        final_max_no_species_plot,
        final_max_total_biomass_plot,
        final_max_no_substrates_plot,
        layout = (2, 2),
        # margin = 5Plots.mm, # Margin around the entire grid
        # spacing = 5Plots.mm, # Spacing between subplots
        plot_title = "final and max values, niche_size$(rl)_$(ru)"
    )
    
        
    max_EUE_times_plot = plot(
        frequencies,
        max_EUE_times,
        xlabel = "Rate", 
        ylabel = "Time",
        legend = false,
        ylims = (0, 3.15e7),
        title = "Time when max EUE was reached, niche_size$(rl)_$(ru)"
    )

    #savefig(final_max_EUEs_plot, joinpath(outdir, "final_max_EUEs_$(frequencies[1])to$(frequencies[end])_rates.png"))
    #savefig(max_EUE_times_plot, joinpath(outdir, "max_EUE_times_$(frequencies[1])to$(frequencies[end])_rates.png"))
    savefig(all_4_final_max_plots, joinpath(outdir, "all_4_final_max_plots_$(frequencies[1])to$(frequencies[end])_rates.png"))
    return (nothing)
end

function EUE_3D_plot()
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

    frequencies = [10, 20, 40, 80, 160, 320, 640]

    println("Compiled and input read in!")
    flush(stdout)

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

        # collect data
        push!(community_EUE_array, community_EUE)

    end

    # Define output directory and if necessary make it
    outdir = joinpath(pwd(), "Output", "niche_size$(rl)_$(ru)", "Immigration_plots")
    mkpath(outdir)
      

    # Initialize the 3D plot with correct setup
    EUE_3D_plot = plot(
        title = "3D Plot of EUE vs Rate vs Time, niche_size$(rl)_$(ru)",
        xlabel = "Time (s)",
        ylabel = "Rate",
        zlabel = "EUE",
        legend = false,
        camera = (60, 30),
        zlims = (0, 1)
    )

    # Add lines for each rate
    for (i, freq) in enumerate(frequencies)
        t_data = t_times
        f_data = fill(freq, length(t_times))
        eue_data = community_EUE_array[i]
        
        plot!(
            EUE_3D_plot,
            t_data, 
            f_data, 
            eue_data, 
            seriestype = :path3d, 
            legend = false,
            camera = (60, 30),
            zlims = (0, 1)
        )
    end

    savefig(EUE_3D_plot, joinpath(outdir, "EUE_3D_plot_$(frequencies[1])to$(frequencies[end])_frequencies.png"))
    return (nothing)
end

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


    num_species_vs_EUE_plot = plot(
        xlabel = "no. species",
        ylabel = "EUE",
        ylims = (0, 1),
        title = "species diversity and EUE, niche_size$(rl)_$(ru)"
    )

    # Add each frequency series to the plot with labels
    for (i, freq) in enumerate(frequencies)   
        plot!(
            num_species_vs_EUE_plot, 
            num_species_array[i], 
            community_EUE_array[i], 
            label="$(freq) events",  
            ms=2, 
            ma=0.5,
            msw = 0,
            seriestype = :scatter
            )
    end

    savefig(num_species_vs_EUE_plot, joinpath(outdir, "num_species_vs_EUE_plot_$(frequencies[1])to$(frequencies[end])_frequencies.png"))

    return (nothing) 
end

@time all_plots()
@time final_max_EUE_plots()
@time EUE_3D_plot()
@time no_species_vs_EUE_plot()

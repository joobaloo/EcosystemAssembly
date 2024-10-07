using JLD2
using Plots
using StatsPlots
using HypothesisTests

"""
    all_plots()

    This function produces 5 time series graphs plotting different variables against time (seconds). 
    The variables are:
        1. Energy use efficiency (EUE), 
        2. Species richness,
        3. Shannon diveristy (H), 
        4. Substrate diversifaction
        5. Total biomass (cells per L)
    
    Inputs:
        This function takes 4 arguments from the command line:
            1. num_immigrants = number of immigration strains per immigration event
            2. rl = lower bound of niche size
            3. ru = upper bound of niche size
            4. rps = number of repeat simulations.

        This function also requires the corresponding RunStats#events_#immigrants.jld files to extract relevant data.
    
    Outputs:
        1. 5 png files for time series plots with each variable
        2. A png file with 4 subplots of each variable (excluding species richness)
        3. A png file comparing species richness and shannon diversity

    Note: 
        Immigration rates are hard-coded as 10, 20, 40, 80, 160 and 320 immigration events per year.
        Simulation length is hardcoded to last a year (3.15e7 seconds).

"""
function all_plots()

        # Preallocate the variables I want to extract from the input
        num_immigrants = 0
        rl = 0
        ru = 0
        rps = 0

        # Check that all arguments can be converted to integers
        try
            num_immigrants = parse(Int64, ARGS[1])
            rl = parse(Int64, ARGS[2])
            ru = parse(Int64, ARGS[3])
            rps = parse(Int64, ARGS[4])
        catch e
            error("Need to provide an integer")
        end
    
        println("Compiled and input read in!")
        flush(stdout)

        # Define immigration rates and simulation length in seconds
        frequencies = [10, 20, 40, 80, 160, 320]
        sim_length = 3.15e7 # 1 year
        
        # Open the JLD file and load the time data while checking it exists
        data_dir = joinpath(pwd(), "Output", "niche_size$(rl)_$(ru)", "$(num_immigrants)immigrants", "$(frequencies[1])events")
        stats_file = joinpath(data_dir, "RunStats$(frequencies[1])events_$(num_immigrants)immigrants.jld")
        if !isfile(stats_file)
            error("missing stats file for $(frequencies[1]) events 1 immigrant simulations")
        end
        t_times = load(stats_file, "times")

        # Initialise variable arrays
        community_EUE_array = []
        num_species_array = []
        num_substrates_array = []
        total_biomass_array = []
        shannon_array =[]

        # Open the JLD file and load variable data for each immigration rate
        for i in frequencies
            data_dir = joinpath(
            pwd(), "Output", "niche_size$(rl)_$(ru)", "$(num_immigrants)immigrants", "$(i)events")
            stats_file = joinpath(data_dir, "RunStats$(i)events_$(num_immigrants)immigrants.jld")
            if ~isfile(stats_file)
                error("missing stats file for $(i)events_$(num_immigrants)immigrants simulations")
            end
        
            # load simulation data
            t_times = load(stats_file, "times")
            community_EUE = load(stats_file, "mean_community_EUE")
            num_species = load(stats_file, "mean_surviving_species")
            num_substrates = load(stats_file, "mean_no_substrates")
            total_biomass = load(stats_file, "mean_total_biomass_of_viable_species")
            shannon = load(stats_file, "mean_shannon_diversity")
            
            # collect data
            push!(community_EUE_array, community_EUE)
            push!(num_species_array, num_species)
            push!(num_substrates_array, num_substrates)
            push!(total_biomass_array, total_biomass)
            push!(shannon_array, shannon)
        end

        # Define output directory and if necessary make it
        outdir = joinpath(pwd(), "Output", "niche_size$(rl)_$(ru)", "Immigration_plots")
        mkpath(outdir)
        
        # Define colour palette
        colour_palette = cgrad(:blues, length(frequencies))

        # Initialise plots
        EUE_plot = plot(
            margin = 10Plots.mm,
            #legend = false,
            xlabel="Time (s)",
            ylabel="EUE",
            tickfontsize = 12,
            )
        num_species_plot = plot(
            #legend = false,
            xlabel="Time (s)",
            tickfontsize = 12,
            ylabel="Number of Species"
            )
        num_substrates_plot = plot(
            #legend = false,
            xlabel="Time",
            tickfontsize = 12,
            ylabel="Substrate Diversification"
            )
        biomass_plot = plot(
            #legend = false,
            xlabel="Time (s)",
            tickfontsize = 12,
            ylabel="Total Biomass (cells per L)"
            )
        shannon_plot = plot(
            #legend = false,
            xlabel="Time (s)",
            tickfontsize = 12,
            ylabel="Shannon Diversity Index (H)"
        )
        
        # Add each immigration rate series to the plot with labels
        for (i, freq) in enumerate(frequencies)
            plot!(
                EUE_plot, 
                t_times, 
                xlims = (0, sim_length),
                ylim = (0,1), 
                community_EUE_array[i],
                tickfontsize = 12,
                color = colour_palette[i], 
                label="$(freq) events"
            )
            plot!(
                num_species_plot, 
                t_times, 
                xlims = (0, sim_length), 
                num_species_array[i],
                ylim = (0, ceil((maximum(num_species_array[i])*1.2))),
                color = colour_palette[i],
                label="$(freq) events"
            )
            plot!(
                num_substrates_plot, 
                t_times, 
                xlims = (0, sim_length), 
                ylim = (0, ceil((maximum(num_substrates_array[i])*1.2))),  
                tickfontsize = 12,
                num_substrates_array[i], 
                color = colour_palette[i],
                label="$(freq) events"
            )
            plot!(
                biomass_plot, 
                t_times, 
                xlims = (0, sim_length),
                ylim = (0, ceil((maximum(total_biomass_array[i])*1.2))),  
                tickfontsize = 12,
                total_biomass_array[i], 
                color = colour_palette[i],
                label="$(freq) events"
            )
            plot!(
                shannon_plot, 
                t_times, 
                xlims = (0, sim_length), 
                shannon_array[i], 
                ylim = (0, ceil((maximum(shannon_array[i])*1.2))), 
                tickfontsize = 12,
                color = colour_palette[i],
                label="$(freq) events"
            )
        end

        # Plot EUE, shannon diversity, substrate diversifaction and total biomass in 1 plot
        p = plot(
            EUE_plot, 
            shannon_plot, 
            num_substrates_plot, 
            biomass_plot, 
            layout = (2, 2), 
            plot_title = "niche_size:$(rl)-$(ru), $(frequencies[1])to$(frequencies[end])_rates, 1 immigrant, $(rps) repeats",
            size = (1200, 800),
            margin = 10Plots.mm
            )
        
        # Plot comparison of species richness and shannon diversity
        compare_diversity = plot(
            num_species_plot,
            shannon_plot,
            layout = (1,2),
            size = (1200, 400),
            margin = 10Plots.mm
            )

        savefig(p, joinpath(outdir, "all_4_plots_$(frequencies[1])to$(frequencies[end])_frequencies.png"))
        savefig(EUE_plot, joinpath(outdir, "EUE_$(frequencies[1])to$(frequencies[end])_frequencies.png"))
        savefig(shannon_plot, joinpath(outdir, "shannon_$(frequencies[1])to$(frequencies[end])_frequencies.png"))
        savefig(num_substrates_plot, joinpath(outdir, "num_substrates_$(frequencies[1])to$(frequencies[end])_frequencies.png"))
        savefig(num_species_plot, joinpath(outdir, "num_species_$(frequencies[1])to$(frequencies[end])_frequencies.png"))
        savefig(biomass_plot, joinpath(outdir, "biomass_$(frequencies[1])to$(frequencies[end])_frequencies.png"))
        savefig(compare_diversity , joinpath(outdir, "shannon_or_richness_$(frequencies[1])to$(frequencies[end])_frequencies.png"))
        return (nothing)
end

"""
    final_max_EUE_plots()

    This function produces 5 plots showing how the final and max values of 5 varaible change as immigration rate increases. 
    The variables are:
        1. Energy use efficiency (EUE), 
        2. Species richness,
        3. Shannon diveristy (H), 
        4. Substrate diversifaction
        5. Total biomass (cells per L)
    
    Inputs:
        This function takes 4 arguments from the command line:
            1. num_immigrants = number of immigration strains per immigration event
            2. rl = lower bound of niche size
            3. ru = upper bound of niche size
            4. rps = number of repeat simulations.

        This function also requires the corresponding RunStats#events_#immigrants.jld files to extract relevant data.
    
    Outputs:
        1. 5 png files for the plots with each variable
        2. A png file with 4 subplots of each variable (excluding species richness)

    Note: 
        Immigration rates are hard-coded as 10, 20, 40, 80, 160 and 320 immigration events per year.

"""
function final_max_EUE_plots()
     # Preallocate the variables I want to extract from the input
     num_immigrants = 0
     rl = 0
     ru = 0
     rps = 0

     # Check that all arguments can be converted to integers
     try
         num_immigrants = parse(Int64, ARGS[1])
         rl = parse(Int64, ARGS[2])
         ru = parse(Int64, ARGS[3])
         rps = parse(Int64, ARGS[4])
     catch e
         error("Need to provide an integer")
     end
 
     println("Compiled and input read in!")
     flush(stdout)

     # Define immigration rates
     frequencies = [10, 20, 40, 80, 160, 320]
     
     # Open the JLD file and load the time data while checking it exists
     data_dir = joinpath(pwd(), "Output", "niche_size$(rl)_$(ru)", "$(num_immigrants)immigrants", "$(frequencies[1])events")
     stats_file = joinpath(data_dir, "RunStats$(frequencies[1])events_$(num_immigrants)immigrants.jld")
     if !isfile(stats_file)
         error("missing stats file for $(frequencies[1]) events 1 immigrant simulations")
     end
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

    final_shannon = []
    final_shannon_SE = []
    max_shannon = []
    max_shannon_SE = []


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
        shannon = load(stats_file, "mean_shannon_diversity")
        shannon_SD = load(stats_file, "sd_shannon_diversity")
    
        max_EUE_value, max_EUE_index = findmax(filter(!isnan, community_EUE))
        max_surviving_species_value, max_surviving_species_index = findmax(filter(!isnan, surviving_species))
        max_total_biomass_value, max_total_biomass_index = findmax(filter(!isnan, total_biomass))
        max_no_substrates_value, max_no_substrates_index = findmax(filter(!isnan, no_substrates))
        max_shannon_value, max_shannon_index = findmax(filter(!isnan, shannon))

        # collect data
        factor = (20-1)/3
        push!(community_EUE_array, community_EUE)
        push!(max_EUE_times, t_times[max_EUE_index])

        push!(final_EUE, community_EUE[end])
        push!(final_EUE_SE, community_EUE_SD[end]/factor)
        push!(max_EUE, max_EUE_value)
        push!(max_EUE_SE, community_EUE_SD[max_EUE_index]/factor)
        
        push!(final_surviving_species, surviving_species[end])
        push!(final_surviving_species_SE, surviving_species_SD[end]/factor)
        push!(max_surviving_species, max_surviving_species_value)
        push!(max_surviving_species_SE, surviving_species_SD[max_surviving_species_index]/factor)

        push!(final_total_biomass, total_biomass[end])
        push!(final_total_biomass_SE, total_biomass_SD[end]/factor)
        push!(max_total_biomass, max_total_biomass_value)
        push!(max_total_biomass_SE, total_biomass_SD[max_total_biomass_index]/factor)

        push!(final_no_substrates, no_substrates[end])
        push!(final_no_substrates_SE, no_substrates_SD[end]/factor)
        push!(max_no_substrates, max_no_substrates_value)
        push!(max_no_substrates_SE, no_substrates_SD[max_no_substrates_index]/factor)

        push!(final_shannon, shannon[end])
        push!(final_shannon_SE, shannon_SD[end]/factor)
        push!(max_shannon, max_shannon_value)
        push!(max_shannon_SE, shannon_SD[max_shannon_index]/factor)
    end

    # Define output directory and if necessary make it
    outdir = joinpath(pwd(), "Output", "niche_size$(rl)_$(ru)", "Immigration_plots")
    mkpath(outdir)

    EUE_matrix = hcat(final_EUE, max_EUE)        
    EUE_SE_matrix = hcat(final_EUE_SE, max_EUE_SE)
    no_species_matrix = hcat(final_surviving_species, max_surviving_species)        
    no_species_SE_matrix = hcat(final_surviving_species_SE, max_surviving_species_SE)
    total_biomass_matrix = hcat(final_total_biomass, max_total_biomass)        
    total_biomass_SE_matrix = hcat(final_total_biomass_SE, max_total_biomass_SE)
    no_substrates_matrix = hcat(final_no_substrates, max_no_substrates)        
    no_substrates_SE_matrix = hcat(final_no_substrates_SE, max_no_substrates_SE)
    shannon_matrix = hcat(final_shannon, max_shannon)        
    shannon_SE_matrix = hcat(final_shannon_SE, max_shannon_SE)
    
    final_max_EUEs_plot = plot(
        EUE_matrix,
        ribbon = EUE_SE_matrix,
        label = ["Final" "Max"],
        xlabel="Rate of Immigration",
        xticks = (1:length(frequencies), frequencies),
        ylabel ="EUE",
        ylims = (0, 1),
        tickfontsize = 12,
        linewidth = 3,
        #legend = false
        #title = "final and max EUE, niche_size$(rl)_$(ru)" 
    )
    
    final_max_no_species_plot = plot(
        no_species_matrix,
        ribbon = no_species_SE_matrix,
        label = ["final" "max"],
        #xlabel="Rates",
        xticks = (1:length(frequencies), frequencies),
        ylabel ="no. species",
        ylim = (0,32),
        #legend = false
        #title = "final and max no. species, niche_size$(rl)_$(ru)" 
    )

    final_max_total_biomass_plot = plot(
        total_biomass_matrix,
        ribbon = total_biomass_SE_matrix,
        label = ["Final" "Max"],
        xlabel=" Immigration rate",
        xticks = (1:length(frequencies), frequencies),
        ylabel="Total biomass (cells per L)",
        ylims = (0, 3.6e14),
        #legend = false
        #title = "final and max total biomass, niche_size$(rl)_$(ru)" 
    )
    
    final_max_no_substrates_plot = plot(
        no_substrates_matrix,
        ribbon = no_substrates_SE_matrix,
        label = ["Final" "Max"],
        xlabel=" Immigration rate",
        xticks = (1:length(frequencies), frequencies),
        ylabel="Substrate Diversification",
        ylims = (0,20),
        #legend = false
        #title = "final and max no. substrates, niche_size$(rl)_$(ru)" 
    )

    final_max_shannon_plot = plot(
        shannon_matrix,
        ribbon = shannon_SE_matrix,
        label = ["Final" "Max"],
        xlabel=" Immigration rate",
        xticks = (1:length(frequencies), frequencies),
        ylabel="Shannon Diversity Index (H)",
        ylims = (0, 2.2),
        #legend = false
        #title = "final and max no. substrates, niche_size$(rl)_$(ru)" 
    )

    all_4_final_max_plots = plot(
        final_max_EUEs_plot,
        final_max_shannon_plot,
        final_max_total_biomass_plot,
        final_max_no_substrates_plot,
        layout = (2, 2),
        size = (1200, 800),
        margin = 10Plots.mm,
        plot_title = "final and max values, niche_size$(rl)_$(ru)"
    )
    
    savefig(final_max_EUEs_plot, joinpath(outdir, "final_max_EUEs_$(frequencies[1])to$(frequencies[end])_rates.png"))
    savefig(final_max_shannon_plot, joinpath(outdir, "final_max_shannon_$(frequencies[1])to$(frequencies[end])_rates.png"))
    savefig(final_max_no_substrates_plot, joinpath(outdir, "final_max_no_substrates_$(frequencies[1])to$(frequencies[end])_rates.png"))
    savefig(final_max_no_species_plot, joinpath(outdir, "final_max_no_species_$(frequencies[1])to$(frequencies[end])_rates.png"))
    savefig(final_max_total_biomass_plot, joinpath(outdir, "final_max_total_biomass_$(frequencies[1])to$(frequencies[end])_rates.png"))
    savefig(all_4_final_max_plots, joinpath(outdir, "all_4_final_max_plots_$(frequencies[1])to$(frequencies[end])_rates.png"))
    return (nothing)
end

"""
    EUE_3D_plot()

    This function creates a 3D plot with energy use efficiency (EUE), time (seconds) and immigration rate.
    
    Inputs:
        This function takes 4 arguments from the command line:
            1. num_immigrants = number of immigration strains per immigration event
            2. rl = lower bound of niche size
            3. ru = upper bound of niche size
            4. rps = number of repeat simulations.

        This function also requires the corresponding RunStats#events_#immigrants.jld files to extract relevant data.
    
    Outputs:
        This function outputs one png file of the 3D plot

    Note: 
        Immigration rates are hard-coded as 10, 20, 40, 80, 160 and 320 immigration events per year.

"""
function EUE_3D_plot()

    # Preallocate the variables I want to extract from the input
    num_immigrants = 0
    rl = 0
    ru = 0
    rps = 0

    # Check that all arguments can be converted to integers
    try
        num_immigrants = parse(Int64, ARGS[1])
        rl = parse(Int64, ARGS[2])
        ru = parse(Int64, ARGS[3])
        rps = parse(Int64, ARGS[4])
    catch e
        error("Need to provide an integer")
    end

    println("Compiled and input read in!")
    flush(stdout)

    # Define immigration rates and simulation length in seconds
    frequencies = [10, 20, 40, 80, 160, 320]
    
    # Open the JLD file and load the time data while checking it exists
    data_dir = joinpath(pwd(), "Output", "niche_size$(rl)_$(ru)", "$(num_immigrants)immigrants", "$(frequencies[1])events")
    stats_file = joinpath(data_dir, "RunStats$(frequencies[1])events_$(num_immigrants)immigrants.jld")
    if !isfile(stats_file)
        error("missing stats file for $(frequencies[1]) events 1 immigrant simulations")
    end
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
        #legend = false,
        camera = (60, 30),
        label = ["10" "20" "40" "80" "160" "320"],
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

#@time all_plots()
#@time final_max_EUE_plots()
@time EUE_3D_plot()


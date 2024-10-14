using JLD2
using Plots

function plot_variation()
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
        shannon_array = []
        num_substrates_array = []
        total_biomass_array = []
        
        community_EUE_se_array = []
        shannon_se_array = []
        num_substrates_se_array = []
        total_biomass_se_array = []

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
            shannon = load(stats_file, "mean_shannon_diversity")
            num_substrates = load(stats_file, "mean_no_substrates")
            total_biomass = load(stats_file, "mean_total_biomass_of_viable_species")

            community_EUE_sd = load(stats_file, "sd_community_EUE")
            shannon_sd = load(stats_file, "sd_shannon_diversity")
            num_substrates_sd = load(stats_file, "sd_no_substrates")
            total_biomass_sd = load(stats_file, "sd_total_biomass_of_viable_species")

            # collect data
            push!(community_EUE_array, community_EUE)
            push!(shannon_array, shannon)
            push!(num_substrates_array, num_substrates)
            push!(total_biomass_array, total_biomass)

            push!(community_EUE_se_array, community_EUE_sd/3)
            push!(shannon_se_array, shannon_sd/3)
            push!(num_substrates_se_array, num_substrates_sd/3)
            push!(total_biomass_se_array, total_biomass_sd/3)

        end

        # Define output directory and if necessary make it
        outdir = joinpath(pwd(), "Output", "niche_size$(rl)_$(ru)", "Immigration_plots", "variation_plots")
        mkpath(outdir)
        custom_layout = @layout [a b; c d; e f ]
        # initialise plots
        EUE_plot = plot(
            legend = false,
            layout = custom_layout,
            size=(1500,2000),
            plot_title = "EUE over time, niche_size$(rl)_$(ru)",
            xlabel="Time", 
            ylabel="EUE",
            margin = 10Plots.mm
            
            )
        shannon_plot = plot(
            legend = false,
            layout = custom_layout,
            size=(2000,1000),
            plot_title = "shannon Diversity over time, niche_size$(rl)_$(ru)",
            xlabel="Time",
            ylabel="Number of species",
            margin = 10Plots.mm
            )
        num_substrates_plot = plot(
            legend = false,
            layout = custom_layout,
            size=(2000,1000),
            plot_title = "Number of substrates over time, niche_size$(rl)_$(ru)",
            xlabel="Time",
            ylabel="Number of substrates",
            margin = 10Plots.mm
            )
        biomass_plot = plot(
            legend = false,
            layout = custom_layout,
            size=(2000,1000),
            plot_title = "Total biomass over time, niche_size$(rl)_$(ru)",
            xlabel="Time",
            ylabel="Total Biomass",
            margin = 10Plots.mm
            )

        # Add each frequency series to the plot with labels
        for (i, freq) in enumerate(frequencies)
            plot!(
                EUE_plot, 
                t_times, 
                community_EUE_array[i],
                ribbon = community_EUE_se_array[i],
                xlims = (0, 6.3e7 / 2),
                ylim = (0,1), 
                subplot = i,
                title="$(freq)/year"
            )
            plot!(
                shannon_plot, 
                t_times, 
                shannon_array[i],
                ribbon = shannon_se_array[i],
                xlims = (0, 6.3e7 / 2),
                ylim = (0,32),
                subplot = i,
                label="$(freq) events"
            )
            plot!(
                num_substrates_plot, 
                t_times, 
                num_substrates_array[i], 
                ribbon = num_substrates_se_array[i],
                xlims = (0, 6.3e7 / 2),
                subplot = i,
                ylims = (0,20),  
                label="$(freq) events"
            )
            plot!(biomass_plot, 
                t_times, 
                total_biomass_array[i],
                ribbon = total_biomass_se_array[i],
                ylims = (0, 2.3e14),
                subplot = i,
                label="$(freq) events"
            )
        end

        savefig(EUE_plot, joinpath(outdir, "EUE_$(frequencies[1])to$(frequencies[end])_frequencies.png"))
        savefig(shannon_plot, joinpath(outdir, "shannon_$(frequencies[1])to$(frequencies[end])_frequencies.png"))
        savefig(num_substrates_plot, joinpath(outdir, "num_substrates_$(frequencies[1])to$(frequencies[end])_frequencies.png"))
        savefig(biomass_plot, joinpath(outdir, "biomass_$(frequencies[1])to$(frequencies[end])_frequencies.png"))
        return (nothing)
end

@time plot_variation()
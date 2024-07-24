using JLD2
using Plots

function all_4_plots()
        # Preallocate the variables I want to extract from the input
        num_immigrants = 0
    
        # Check that all arguments can be converted to integers
        try
            num_immigrants = parse(Int64, ARGS[1])
        catch e
            error("Need to provide an integer")
        end
    
        println("Compiled and input read in!")
        flush(stdout)

        # Open the JLD file and load the time data
        data_dir = joinpath(pwd(), "Output", "1immigrants", "50events")

        stats_file = joinpath(data_dir, "RunStats50events_1immigrants.jld")
        
        # Check it actually exists
        if !isfile(stats_file)
            error("missing stats file for 50 events 1 immigrant simulations")
        end

        # Extract time data
        t_times = load(stats_file, "times")

        frequencies = [50, 75, 100, 125, 150]
        community_EUE_array = []
        num_species_array = []
        num_substrates_array = []
        total_biomass_array = []

        for i in frequencies
            # Open the JLD file and load the surviving species data
            data_dir = joinpath(
            pwd(), "Output", "$(num_immigrants)immigrants", "$(i)events")
        
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
        outdir = joinpath(pwd(), "Output", "Immigration_plots")
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
            plot_title = "50 - 150 events, 1 immigrant, 10 repeats",
            size = (1200, 800),
            margin = 10Plots.mm
            )

        savefig(p, joinpath(outdir, "all_4_plots_50to150_frequencies.png"))
        savefig(EUE_plot, joinpath(outdir, "EUE_50to150_frequencies.png"))
        savefig(num_species_plot, joinpath(outdir, "num_species_50to150_frequencies.png"))
        savefig(num_substrates_plot, joinpath(outdir, "num_substrates_50to150_frequencies.png"))
        savefig(biomass_plot, joinpath(outdir, "biomass_50to150_frequencies.png"))
        return (nothing)
end

@time all_4_plots()
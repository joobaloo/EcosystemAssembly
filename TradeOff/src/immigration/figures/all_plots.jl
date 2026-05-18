using JLD2
using Plots
using StatsPlots
using HypothesisTests
using Printf

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
            1. rps = number of repeat simulations.
            2. num_immigrants = number of immigration strains per immigration event
            3. rl = lower bound of niche size
            4. ru = upper bound of niche size

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
        rps = 0
        num_immigrants = 0
        rl = 0
        ru = 0

        # Check that all arguments can be converted to integers
        try
            rps = parse(Int64, ARGS[1])
            num_immigrants = parse(Int64, ARGS[2])
            rl = parse(Int64, ARGS[3])
            ru = parse(Int64, ARGS[4])
        catch e
            error("Need to provide an integer")
        end
    
        println("Compiled and input read in!")
        flush(stdout)

        # Define immigration rates and simulation length in seconds
        #frequencies = [10, 20, 40, 80, 160, 320, 640]
        frequencies = [10, 80, 320]
        sim_length = 3.15e7 * 32 # 4 years

        # Initialise variable arrays
        community_EUE_array = []
        num_species_array = []
        num_substrates_array = []
        total_biomass_array = []
        shannon_array =[]
        t_times_array = []

        # Initialise variable arrays for standard deviations
        community_EUE_sd_array = []
        num_species_sd_array = []
        num_substrates_sd_array = []
        total_biomass_sd_array = []
        shannon_sd_array = []

        # Open the JLD file and load variable data for each immigration rate
        for i in frequencies
            data_dir = joinpath(
            pwd(), "Output", "niche_size$(rl)_$(ru)", "$(num_immigrants)immigrants", "$(i)_a_year_rate")
            stats_file = joinpath(data_dir, "RunStats$(i)_a_year_rate_$(num_immigrants)immigrants.jld2")
            if ~isfile(stats_file)
                error("missing stats file for $(i)events_$(num_immigrants)immigrants simulations")
            end
        
            # Load mean data
            push!(community_EUE_array, load(stats_file, "mean_community_EUE"))
            push!(num_species_array, load(stats_file, "mean_surviving_species"))
            push!(num_substrates_array, load(stats_file, "mean_no_substrates"))
            push!(total_biomass_array, load(stats_file, "mean_total_biomass_of_viable_species"))
            push!(shannon_array, load(stats_file, "mean_shannon_diversity"))
            push!(t_times_array, load(stats_file, "times"))

            # Load SD data
            push!(community_EUE_sd_array, load(stats_file, "sd_community_EUE"))
            push!(num_species_sd_array, load(stats_file, "sd_surviving_species"))
            push!(num_substrates_sd_array, load(stats_file, "sd_no_substrates"))
            push!(total_biomass_sd_array, load(stats_file, "sd_total_biomass_of_viable_species"))
            push!(shannon_sd_array, load(stats_file, "sd_shannon_diversity"))
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
            xlabel="Time (years)",
            xlims=(0,1),
            xticks=0:0.5:1,
            ylabel="EUE",
            tickfontsize = 12,
            )
        num_species_plot = plot(
            #legend = false,
            xlabel="Time (years)",
            xlims=(0,1),
            xticks=0:0.5:1,
            tickfontsize = 12,
            ylabel="Number of Species"
            )
        num_substrates_plot = plot(
            #legend = false,
            xlabel="Time (years)",
            xlims=(0,1),
            xticks=0:0.5:1,
            tickfontsize = 12,
            ylabel="Substrate Diversification"
            )
        biomass_plot = plot(
            #legend = false,
            xlabel="Time (years)",
            xlims=(0,1),
            xticks=0:0.5:1,
            tickfontsize = 12,
            ylabel="Total Biomass (cells per L)"
            )
        shannon_plot = plot(
            #legend = false,
            xlabel="Time (years)",
            xlims=(0,1),
            xticks=0:0.5:1,
            tickfontsize = 12,
            ylabel="Shannon Diversity Index (H)"
        )

        clean_sd(sd) = replace(sd, NaN => 0.0)
        
        # Add each immigration rate series to the plot with labels
        for (i, freq) in enumerate(frequencies)
            # convert seconds → years
            scaled_time = t_times_array[i] ./ 3.15e7 
            plot!(
                EUE_plot, 
                scaled_time, 
                ylim = (0,1), 
                community_EUE_array[i],
                ribbon = clean_sd(community_EUE_sd_array[i]), fillalpha = 0.2,
                tickfontsize = 12,
                color = colour_palette[i], 
                label="$(freq) events"
            )
            plot!(
                num_species_plot, 
                scaled_time, 
                num_species_array[i],
                ribbon = clean_sd(num_species_sd_array[i]), fillalpha = 0.2,
                ylim = (0, ceil((maximum(num_species_array[i])*1.2))),
                color = colour_palette[i],
                label="$(freq) events"
            )
            plot!(
                num_substrates_plot, 
                scaled_time,
                ylim = (0, ceil((maximum(num_substrates_array[i])*1.2))),  
                tickfontsize = 12,
                num_substrates_array[i], 
                ribbon = clean_sd(num_substrates_sd_array[i]), fillalpha = 0.2,
                color = colour_palette[i],
                label="$(freq) events"
            )
            plot!(
                biomass_plot, 
                scaled_time, 
                ylim = (0, ceil((maximum(total_biomass_array[i])*1.2))),  
                tickfontsize = 12,
                total_biomass_array[i], 
                ribbon = clean_sd(total_biomass_sd_array[i]), fillalpha = 0.2,
                color = colour_palette[i],
                label="$(freq) events"
            )
            plot!(
                shannon_plot, 
                scaled_time, 
                shannon_array[i], 
                ribbon = clean_sd(shannon_sd_array[i]), fillalpha = 0.2,
                ylim = (0, ceil((maximum(shannon_array[i])*1.2))), 
                tickfontsize = 12,
                color = colour_palette[i],
                label="$(freq) events"
            )
        end

        # Plot EUE, species richness, substrate diversifaction and total biomass in 1 plot
        p = plot(
            EUE_plot, 
            num_species_plot, 
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
            size = (2000, 1200),
            margin = 10Plots.mm
            )

        savefig(p, joinpath(outdir, "all_4_plots_$(frequencies[1])to$(frequencies[end])_frequencies.png"))
        # savefig(EUE_plot, joinpath(outdir, "EUE_$(frequencies[1])to$(frequencies[end])_frequencies.png"))
        # savefig(shannon_plot, joinpath(outdir, "shannon_$(frequencies[1])to$(frequencies[end])_frequencies.png"))
        # savefig(num_substrates_plot, joinpath(outdir, "num_substrates_$(frequencies[1])to$(frequencies[end])_frequencies.png"))
        # savefig(num_species_plot, joinpath(outdir, "num_species_$(frequencies[1])to$(frequencies[end])_frequencies.png"))
        # savefig(biomass_plot, joinpath(outdir, "biomass_$(frequencies[1])to$(frequencies[end])_frequencies.png"))
        # savefig(compare_diversity , joinpath(outdir, "shannon_or_richness_$(frequencies[1])to$(frequencies[end])_frequencies.png"))
        return (nothing)
end

function all_4_plots_with_SD()

        # Preallocate the variables I want to extract from the input
        rps = 0
        num_immigrants = 0
        rl = 0
        ru = 0

        # Check that all arguments can be converted to integers
        try
            rps = parse(Int64, ARGS[1])
            num_immigrants = parse(Int64, ARGS[2])
            rl = parse(Int64, ARGS[3])
            ru = parse(Int64, ARGS[4])
        catch e
            error("Need to provide an integer")
        end
    
        println("Compiled and input read in!")
        flush(stdout)

        # Define immigration rates and simulation length in seconds
        frequencies = [10, 80, 320]
        sim_length = 3.15e7 * 32 # 4 years

        # Initialise variable arrays for means
        community_EUE_array = []
        num_species_array = []
        num_substrates_array = []
        total_biomass_array = []
        shannon_array = []
        t_times_array = []

        # Initialise variable arrays for standard deviations
        community_EUE_sd_array = []
        num_species_sd_array = []
        num_substrates_sd_array = []
        total_biomass_sd_array = []
        shannon_sd_array = []

        # Open the JLD2 file and load variable data for each immigration rate
        # Note: averages.jl saves as .jld2, ensuring compatibility here
        for i in frequencies
            data_dir = joinpath(
                pwd(), "Output", "niche_size$(rl)_$(ru)", "$(num_immigrants)immigrants", "$(i)_a_year_rate")
            
            # Note: Changed to .jld2 to match the output of averages.jl
            stats_file = joinpath(data_dir, "RunStats$(i)_a_year_rate_$(num_immigrants)immigrants.jld2")
            
            if ~isfile(stats_file)
                error("missing stats file for $(i)events_$(num_immigrants)immigrants simulations")
            end
        
            # Load mean data
            push!(community_EUE_array, load(stats_file, "mean_community_EUE"))
            push!(num_species_array, load(stats_file, "mean_surviving_species"))
            push!(num_substrates_array, load(stats_file, "mean_no_substrates"))
            push!(total_biomass_array, load(stats_file, "mean_total_biomass_of_viable_species"))
            push!(shannon_array, load(stats_file, "mean_shannon_diversity"))
            push!(t_times_array, load(stats_file, "times"))

            # Load SD data
            push!(community_EUE_sd_array, load(stats_file, "sd_community_EUE"))
            push!(num_species_sd_array, load(stats_file, "sd_surviving_species"))
            push!(num_substrates_sd_array, load(stats_file, "sd_no_substrates"))
            push!(total_biomass_sd_array, load(stats_file, "sd_total_biomass_of_viable_species"))
            push!(shannon_sd_array, load(stats_file, "sd_shannon_diversity"))
        end

        # Define output directory
        outdir = joinpath(pwd(), "Output", "niche_size$(rl)_$(ru)", "Immigration_plots")
        mkpath(outdir)
        
        # Define colour palette
        colour_palette = cgrad(:blues, length(frequencies))

        # Initialise plots
        EUE_plot = plot(margin = 10Plots.mm, xlabel="Time (years)", xlims=(0,1), xticks=0:0.5:1, ylabel="EUE", tickfontsize = 12)
        num_species_plot = plot(xlabel="Time (years)", xlims=(0,1), xticks=0:0.5:1, tickfontsize = 12, ylabel="Number of Species")
        num_substrates_plot = plot(xlabel="Time (years)", xlims=(0,1), xticks=0:0.5:1, tickfontsize = 12, ylabel="Substrate Diversification")
        biomass_plot = plot(xlabel="Time (years)", xlims=(0,1), xticks=0:0.5:1, tickfontsize = 12, ylabel="Total Biomass (cells per L)")
        shannon_plot = plot(xlabel="Time (years)", xlims=(0,1), xticks=0:0.5:1, tickfontsize = 12, ylabel="Shannon Diversity Index (H)")
        
        # Add each immigration rate series to the plot with SD ribbons
        for (i, freq) in enumerate(frequencies)
            scaled_time = t_times_array[i] ./ 3.15e7
            
            # Helper to handle NaN/Missing in ribbons if any
            clean_sd(sd) = replace(sd, NaN => 0.0)

            plot!(EUE_plot, scaled_time, community_EUE_array[i], 
                ribbon = clean_sd(community_EUE_sd_array[i]), fillalpha = 0.2,
                color = colour_palette[i], label="$(freq) events", ylim = (0,1))

            plot!(num_species_plot, scaled_time, num_species_array[i], 
                ribbon = clean_sd(num_species_sd_array[i]), fillalpha = 0.2,
                color = colour_palette[i], label="$(freq) events")

            plot!(num_substrates_plot, scaled_time, num_substrates_array[i], 
                ribbon = clean_sd(num_substrates_sd_array[i]), fillalpha = 0.2,
                color = colour_palette[i], label="$(freq) events")

            plot!(biomass_plot, scaled_time, total_biomass_array[i], 
                ribbon = clean_sd(total_biomass_sd_array[i]), fillalpha = 0.2,
                color = colour_palette[i], label="$(freq) events")

            plot!(shannon_plot, scaled_time, shannon_array[i], 
                ribbon = clean_sd(shannon_sd_array[i]), fillalpha = 0.2,
                color = colour_palette[i], label="$(freq) events")
        end

        # Layout and save
        p = plot(EUE_plot, num_species_plot, num_substrates_plot, biomass_plot, 
            layout = (2, 2), size = (1200, 800), margin = 10Plots.mm,
            plot_title = "niche_size:$(rl)-$(ru), $(frequencies[1])to$(frequencies[end])_rates, $(num_immigrants) immigrant, $(rps) repeats")
        
        savefig(p, joinpath(outdir, "all_4_plots_with_SD.png"))
        return (nothing)
end

function all_4_plots_with_SD_3by4()

        # Preallocate variables
        rps = 0
        num_immigrants = 0
        rl = 0
        ru = 0

        # Parse command line arguments
        try
            rps = parse(Int64, ARGS[1])
            num_immigrants = parse(Int64, ARGS[2])
            rl = parse(Int64, ARGS[3])
            ru = parse(Int64, ARGS[4])
        catch e
            error("Need to provide an integer")
        end
    
        println("Compiled and input read in!")
        flush(stdout)

        # Define immigration rates
        frequencies = [10, 80, 320] # [cite: 31]

        # Initialise arrays for data and SDs
        community_EUE_array = []
        num_substrates_array = []
        total_biomass_array = []
        shannon_array = []
        t_times_array = []

        community_EUE_sd_array = []
        num_substrates_sd_array = []
        total_biomass_sd_array = []
        shannon_sd_array = []

        # Load data for each immigration rate
        for i in frequencies
            data_dir = joinpath(pwd(), "Output", "niche_size$(rl)_$(ru)", "$(num_immigrants)immigrants", "$(i)_a_year_rate")
            stats_file = joinpath(data_dir, "RunStats$(i)_a_year_rate_$(num_immigrants)immigrants.jld2") 
            
            if !isfile(stats_file)
                error("missing stats file for $(i)events simulations")
            end
        
            push!(community_EUE_array, load(stats_file, "mean_community_EUE")) 
            push!(num_substrates_array, load(stats_file, "mean_no_substrates")) 
            push!(total_biomass_array, load(stats_file, "mean_total_biomass_of_viable_species")) 
            push!(shannon_array, load(stats_file, "mean_shannon_diversity")) 
            push!(t_times_array, load(stats_file, "times")) 

            push!(community_EUE_sd_array, load(stats_file, "sd_community_EUE")) 
            push!(num_substrates_sd_array, load(stats_file, "sd_no_substrates")) 
            push!(total_biomass_sd_array, load(stats_file, "sd_total_biomass_of_viable_species")) 
            push!(shannon_sd_array, load(stats_file, "sd_shannon_diversity")) 
        end

        outdir = joinpath(pwd(), "Output", "niche_size$(rl)_$(ru)", "Immigration_plots")
        mkpath(outdir)
        
        colour_palette = cgrad(:blues, length(frequencies))
        clean_sd(sd) = replace(sd, NaN => 0.0)

        # Define metrics for the 4 rows
        metric_names = ["EUE", "Substrate Diversification", "Total Biomass (cells/L)", "Shannon Index (H)"]
        data_sets = [community_EUE_array, num_substrates_array, total_biomass_array, shannon_array]
        sd_sets = [community_EUE_sd_array, num_substrates_sd_array, total_biomass_sd_array, shannon_sd_array]

        plot_list = []

        # Build 4x3 grid: Rows = Metrics, Columns = Frequencies
        for m_idx in 1:4
            for f_idx in 1:3
                raw_time = t_times_array[f_idx] # Using raw time (seconds)
                y_data = data_sets[m_idx][f_idx]
                y_sd = sd_sets[m_idx][f_idx]
                
                # Determine plot limits based on raw simulation length
                max_t = maximum(raw_time)

                p = plot(raw_time, y_data, 
                    ribbon = clean_sd(y_sd), 
                    fillalpha = 0.2,
                    color = colour_palette[f_idx],
                    xlims = (0, max_t),
                    ylabel = (f_idx == 1 ? metric_names[m_idx] : ""), # Label y-axis only on leftmost column
                    xlabel = (m_idx == 4 ? "Time (seconds)" : ""),   # Label x-axis only on bottom row
                    title = (m_idx == 1 ? "$(frequencies[f_idx]) events" : ""), # Title only top row
                    legend = false,
                    tickfontsize = 9,
                    guidefontsize = 11
                )
                
                push!(plot_list, p)
            end
        end

        # Arrange in a 4x3 layout
        final_plot = plot(plot_list..., 
            layout = (4, 3), 
            size = (1400, 1200), 
            margin = 10Plots.mm,
            plot_title = "niche_size: $(rl)-$(ru), $(num_immigrants) immigrant, $(rps) repeats")
        
        savefig(final_plot, joinpath(outdir, "multiplot_4x3_raw_time.png"))
        return nothing
end


function all_4_plots_fraction()

        # Preallocate the variables I want to extract from the input
        rps = 0
        num_immigrants = 0
        rl = 0
        ru = 0

        # Check that all arguments can be converted to integers
        try
            rps = parse(Int64, ARGS[1])
            num_immigrants = parse(Int64, ARGS[2])
            rl = parse(Int64, ARGS[3])
            ru = parse(Int64, ARGS[4])
        catch e
            error("Need to provide an integer")
        end
    
        println("Compiled and input read in!")
        flush(stdout)

        # Define immigration rates and simulation length in seconds
        frequencies = [10, 80, 320]
        sim_length = 3.15e7 * 32 # Define the total simulation length 

        # Initialise variable arrays for means
        community_EUE_array = []
        num_species_array = []
        num_substrates_array = []
        total_biomass_array = []
        shannon_array = []
        t_times_array = []

        # Initialise variable arrays for standard deviations
        community_EUE_sd_array = []
        num_species_sd_array = []
        num_substrates_sd_array = []
        total_biomass_sd_array = []
        shannon_sd_array = []

        # Open the JLD2 file and load variable data for each immigration rate
        for i in frequencies
            data_dir = joinpath(
                pwd(), "Output", "niche_size$(rl)_$(ru)", "$(num_immigrants)immigrants", "$(i)_a_year_rate")
            
            stats_file = joinpath(data_dir, "RunStats$(i)_a_year_rate_$(num_immigrants)immigrants.jld2")
            
            if ~isfile(stats_file)
                error("missing stats file for $(i)events_$(num_immigrants)immigrants simulations")
            end
        
            # Load mean and SD data [cite: 35, 36]
            push!(community_EUE_array, load(stats_file, "mean_community_EUE"))
            push!(num_species_array, load(stats_file, "mean_surviving_species"))
            push!(num_substrates_array, load(stats_file, "mean_no_substrates"))
            push!(total_biomass_array, load(stats_file, "mean_total_biomass_of_viable_species"))
            push!(shannon_array, load(stats_file, "mean_shannon_diversity"))
            push!(t_times_array, load(stats_file, "times"))

            push!(community_EUE_sd_array, load(stats_file, "sd_community_EUE"))
            push!(num_species_sd_array, load(stats_file, "sd_surviving_species"))
            push!(num_substrates_sd_array, load(stats_file, "sd_no_substrates"))
            push!(total_biomass_sd_array, load(stats_file, "sd_total_biomass_of_viable_species"))
            push!(shannon_sd_array, load(stats_file, "sd_shannon_diversity"))
        end

        # Define output directory
        outdir = joinpath(pwd(), "Output", "niche_size$(rl)_$(ru)", "Immigration_plots")
        mkpath(outdir)
        
        # Define colour palette
       #colour_palette = cgrad(:blues, length(frequencies))
        colour_palette = palette(:tab10)[1:length(frequencies)]

        # Initialise plots with the new xlabel [cite: 37, 38]
        EUE_plot = plot(margin = 10Plots.mm, xlabel="Time (fraction of simulation)", xlims=(0,1), xticks=0:0.5:1, ylabel="EUE", tickfontsize = 12)
        num_species_plot = plot(xlabel="Time (fraction of simulation)", xlims=(0,1), xticks=0:0.5:1, tickfontsize = 12, ylabel="Number of Species")
        num_substrates_plot = plot(xlabel="Time (fraction of simulation)", xlims=(0,1), xticks=0:0.5:1, tickfontsize = 12, ylabel="Substrate Diversification")
        biomass_plot = plot(xlabel="Time (fraction of simulation)", xlims=(0,1), xticks=0:0.5:1, tickfontsize = 12, ylabel="Total Biomass (cells per L)")
        shannon_plot = plot(xlabel="Time (fraction of simulation)", xlims=(0,1), xticks=0:0.5:1, tickfontsize = 12, ylabel="Shannon Diversity Index (H)")
        
        # Add each immigration rate series to the plot with SD ribbons
        for (i, freq) in enumerate(frequencies)
            # CHANGE: Calculate time as a fraction of total simulation length [cite: 39]
            scaled_time = t_times_array[i] ./ maximum(t_times_array[i])
            
            clean_sd(sd) = replace(sd, NaN => 0.0)

            plot!(EUE_plot, scaled_time, community_EUE_array[i], 
                ribbon = clean_sd(community_EUE_sd_array[i]), fillalpha = 0.2,
                color = colour_palette[i], label="$(freq) events", ylim = (0,1))

            plot!(num_species_plot, scaled_time, num_species_array[i], 
                ribbon = clean_sd(num_species_sd_array[i]), fillalpha = 0.2,
                color = colour_palette[i], label="$(freq) events")

            plot!(num_substrates_plot, scaled_time, num_substrates_array[i], 
                ribbon = clean_sd(num_substrates_sd_array[i]), fillalpha = 0.2,
                color = colour_palette[i], label="$(freq) events")

            plot!(biomass_plot, scaled_time, total_biomass_array[i], 
                ribbon = clean_sd(total_biomass_sd_array[i]), fillalpha = 0.2,
                color = colour_palette[i], label="$(freq) events")

            plot!(shannon_plot, scaled_time, shannon_array[i], 
                ribbon = clean_sd(shannon_sd_array[i]), fillalpha = 0.2,
                color = colour_palette[i], label="$(freq) events")
        end

        # Layout and save [cite: 43]
        p = plot(EUE_plot, num_species_plot, num_substrates_plot, biomass_plot, 
            layout = (2, 2), size = (1200, 800), margin = 10Plots.mm,
            plot_title = "niche_size:$(rl)-$(ru), $(frequencies[1])to$(frequencies[end])_rates, $(num_immigrants) immigrant, $(rps) repeats")
        
        savefig(p, joinpath(outdir, "all_4_plots_fraction.png"))
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
     rps = 0
     num_immigrants = 0
     rl = 0
     ru = 0

     # Check that all arguments can be converted to integers
     try
         rps = parse(Int64, ARGS[1])
         num_immigrants = parse(Int64, ARGS[2])
         rl = parse(Int64, ARGS[3])
         ru = parse(Int64, ARGS[4])
     catch e
         error("Need to provide an integer")
     end
 
     println("Compiled and input read in!")
     flush(stdout)

     # Define immigration rates
     #frequencies = [10, 20, 40, 80, 160, 320, 640]
     frequencies = [10, 80, 320]
     
     # Open the JLD file and load the time data while checking it exists
     data_dir = joinpath(pwd(), "Output", "niche_size$(rl)_$(ru)", "$(num_immigrants)immigrants", "$(frequencies[1])_a_year_rate")
     stats_file = joinpath(data_dir, "RunStats$(frequencies[1])_a_year_rate_$(num_immigrants)immigrants.jld2")
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
        pwd(), "Output", "niche_size$(rl)_$(ru)", "$(num_immigrants)immigrants", "$(i)_a_year_rate")
    
        stats_file = joinpath(data_dir, "RunStats$(i)_a_year_rate_$(num_immigrants)immigrants.jld2")

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
        ylims = (0.4, 0.8),
        tickfontsize = 12,
        linewidth = 3,
        #legend = false
        title = "final and max EUE, niche_size$(rl)_$(ru)" 
    )
    
    final_max_no_species_plot = plot(
        no_species_matrix,
        ribbon = no_species_SE_matrix,
        label = ["final" "max"],
        #xlabel="Rates",
        xticks = (1:length(frequencies), frequencies),
        ylabel ="no. species",
        #ylim = (0,32),
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
        #ylims = (0, 3.6e14),
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
        #ylims = (0,20),
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
        #ylims = (0, 2.2),
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


#
@time all_plots()
#@time all_4_plots_with_SD()
#@time final_max_EUE_plots()
@time all_4_plots_fraction()

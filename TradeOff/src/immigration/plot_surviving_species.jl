using JLD2
using Plots

function plot_surviving_species()

    # Preallocate the variables I want to extract from the input
    rps = 0
    num_immigrations = 0
    num_immigrants = 0

    # Check that all arguments can be converted to integers
    try
        rps = parse(Int64, ARGS[1])
        num_immigrations = parse(Int64, ARGS[2])
        num_immigrants = parse(Int64, ARGS[3])
    catch e
        error("Need to provide 2 integers")
    end

    println("Compiled and input read in!")
    flush(stdout)

    # Open the JLD file and load the surviving species data
    data_dir = joinpath(
    pwd(), "Output", "$(num_immigrants)immigrants", "$(num_immigrations)events")

    stats_file = joinpath(data_dir, "RunStats$(num_immigrations)events_$(num_immigrants)immigrants.jld")
    
    # Check it actually exists
    if ~isfile(stats_file)
        error("missing stats file for $(num_immigrations)events_$(num_immigrants)immigrants simulations")
    end

    # Extract time and surviving species data
    t_times = load(stats_file, "times")
    mean_surviving_species = load(stats_file, "mean_surviving_species")
    sd_surviving_species = load(stats_file, "sd_surviving_species")

    #calculate standard error
    se_surviving_species = sd_surviving_species ./ sqrt.(rps - 1)
    # Define output directory and if necessary make it
    outdir = joinpath(pwd(), "Output", "$(num_immigrants)immigrants", "$(num_immigrations)events")
    mkpath(outdir) 
    
    #println(t_times)
    # Plotting the data
    p = plot(
        t_times,
        mean_surviving_species,
        ribbon = se_surviving_species,
        fillcolor = :lightgreen,
        fillalpha = 0.5, 
        xlabel="Time",
        xlims = (0, 6.3e7*2),
        ylabel="Surviving Species",
        title = "$(num_immigrations) Events, $(num_immigrants) Immigrants, $(rps) repeats",
        legend = false,
        #ylims = (0, 200), 
        # label="Surviving Species over Time", 
        # legend=:topleft
        )
    savefig(p, joinpath(outdir, "$(num_immigrations)events$(num_immigrants)immigrants_surviving_species.png"))
    return (nothing)
end


function plot_time_and_freq()
    # Preallocate the variables I want to extract from the input
    rps = 0
    num_immigrants = 0

    # Check that all arguments can be converted to integers
    try
        rps = parse(Int64, ARGS[1])
        num_immigrants = parse(Int64, ARGS[2])
    catch e
        error("Need to provide 2 integers")
    end

    println("Compiled and input read in!")
    flush(stdout)

    # Open the JLD file and load the time data
    data_dir = joinpath(pwd(), "Output", "1immigrants", "5events")

    stats_file = joinpath(data_dir, "RunStats5events_1immigrants.jld")
    
    # Check it actually exists
    if !isfile(stats_file)
        error("missing stats file for 5 events 1 immigrant simulations")
    end

    # Extract time data
    t_times = load(stats_file, "times")

    frequencies = [10, 20, 30, 40, 50, 60, 70]
    diff_freqs_surviving_species = []

    for i in frequencies
        # Open the JLD file and load the surviving species data
        simulation_dir = joinpath(pwd(), "Output", "$(num_immigrants)immigrants", "$(i)events")
        simulation_stats = joinpath(simulation_dir, "RunStats$(i)events_$(num_immigrants)immigrants.jld")

        # Check it actually exists
        if !isfile(simulation_stats)
            error("missing stats file for $(i) events $(num_immigrants) immigrants simulations")
        end

        # Load mean surviving species
        column = load(simulation_stats, "mean_surviving_species")
        
        # Collect the data
        push!(diff_freqs_surviving_species, column)
    end

    # Define output directory and if necessary make it
    outdir = joinpath(pwd(), "Output", "Immigration_plots")
    mkpath(outdir)

    # Plotting the data
    p = plot(xlabel="Time", ylabel="Surviving Species", title="5 - 35 Events, $(num_immigrants) Immigrants, $(rps) repeats")

    # Add each frequency series to the plot with labels
    for (i, freq) in enumerate(frequencies)
        plot!(p, t_times, xlims = (0, 6.3e7 * 2), diff_freqs_surviving_species[i], label="$(freq) events")
    end

    savefig(p, joinpath(outdir, "10to70_frequencies_$(num_immigrants)immigrants_surviving_species.png"))
    return nothing
end



@time plot_surviving_species()
#@time plot_time_and_freq()

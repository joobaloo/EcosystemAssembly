using JLD
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
        pwd(), "Output", "$(num_immigrations)events_$(num_immigrants)immigrants")

    stats_file = joinpath(data_dir, "RunStats$(num_immigrations)events_$(num_immigrants)immigrants.jld")
    
    # Check it actually exists
    if ~isfile(stats_file)
        error("missing stats file for $(num_immigrations)events_$(num_immigrants)immigrants simulations")
    end

    # Extract time and surviving species data
    t_times = load(stats_file, "times")
    mean_surviving_species = load(stats_file, "mean_surviving_species")
    #println(typeof(mean_surviving_species))
    #println(mean_surviving_species)
    # Define output directory and if necessary make it
    outdir = joinpath(pwd(), "Output", "Imm_plots")
    mkpath(outdir) 
    
    #println(t_times)
    # Plotting the data
    p = plot(
        t_times,
        mean_surviving_species, 
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

@time plot_surviving_species()
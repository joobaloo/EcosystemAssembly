using JLD2
using Plots

function plot_species_EUE()

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

    stats_file = joinpath(data_dir,  "AvRun$(rps)Data.jld")
    
    # Check it actually exists
    if ~isfile(stats_file)
        error("missing stats file for $(num_immigrations)events_$(num_immigrants)immigrants simulations")
    end

    # Extract time and species EUE data
    t_times = load(stats_file, "T")
    species_EUEs = load(stats_file, "species_EUEs")


    # Define output directory and if necessary make it
    outdir = joinpath(pwd(), "Output", "$(num_immigrants)immigrants", "$(num_immigrations)events")
    mkpath(outdir) 
    
    # Plotting the data
    p = plot(
        t_times,
        species_EUEs,
        xlabel="Time",
        xlims = (0, 6.3e7*2),
        ylabel="EUE",
        title = "Species EUE $(num_immigrations) Events, $(num_immigrants) Immigrants, $(rps) repeats",
        #legend = false,
        #ylims = (0, 200), 
        # label="Surviving Species over Time", 
        legend= false
        )
    savefig(p, joinpath(outdir, "$(num_immigrations)events$(num_immigrants)immigrants_species_EUEs.png"))
    return (nothing)
end

function plot_community_EUE()
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
            error("Need to provide 3 integers")
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
        

        # Extract time and community EUE data
        t_times = load(stats_file, "times")
        community_EUE = load(stats_file, "mean_community_EUE")
        
         # Verify that the lengths of t_times and community_EUE are the same
        if length(t_times) != length(community_EUE)
            error("Length mismatch: t_times has length $(length(t_times)), but community_EUE has length $(length(community_EUE)).")
        end

        # Define output directory and if necessary make it
        outdir = joinpath(pwd(), "Output", "$(num_immigrants)immigrants", "$(num_immigrations)events")
        mkpath(outdir) 
        
        # Plotting the data
        p = plot(
            t_times,
            community_EUE,
            xlabel="Time",
            #xlims = (0, 6.3e7*2),
            ylabel="EUE",
            ylims = (0, 1),
            title = "Communtiy EUE $(num_immigrations) Events, $(num_immigrants) Immigrants, $(rps) repeats",
            legend = false,
            )
        savefig(p, joinpath(outdir, "$(num_immigrations)events$(num_immigrants)immigrants_community_EUE.png"))
        return (nothing)
end


#time plot_species_EUE()
@time plot_community_EUE()
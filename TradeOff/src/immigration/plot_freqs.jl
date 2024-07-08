using TradeOff
using Plots
using JLD

include("../immigration/simulation_functions.jl")

function frequency_plots()
    # Preallocate the variables I want to extract from the input
    rps = 0
    sim_type = 0
    total_time = 6.3e7 # this is approx 2 years in seconds
    num_immigrations = 0
    num_immigrants = 0

    # Check that all arguments can be converted to integers
    try
        rps = parse(Int64, ARGS[1])
        sim_type = parse(Int64, ARGS[2])
        num_immigrations = parse(Int64, ARGS[3])
        num_immigrants = parse(Int64, ARGS[4])
    catch e
        error("Need to provide 4 integers")
    end

    if sim_type != 1
        error("Only use simulation type 1")
    end

    println("Compiled and input read in!")
    flush(stdout)
    
    # Extract other simulation parameters from the fuLoadError: MethodError: no method matching merge_datanction
    #Np, Nt, M, d, μrange = imm_sim_paras(sim_type)
    #println(μrange)

    # Read in appropriate files
    pfile = "Output/$(num_immigrations)events_$(num_immigrants)immigrants/Parameters.jld"
    #println("Output/$(num_immigrations)events_$(num_immigrants)immigrants/Parameters.jld")
    if ~isfile(pfile)
        error("$(num_immigrants) immigrations run $(rps) is missing a parameter file")
    end
    
    ofile = "Output/$(num_immigrations)events_$(num_immigrants)immigrants/Run$(rps)Data.jld"
    if ~isfile(ofile)
        error("$(num_immigrations) immigration events with (num_immigrants)immigrants run $(rps) is missing an output file")
    end

    # Read in relevant data
    ps = load(pfile, "ps")
    traj = load(ofile, "traj")
    T = load(ofile, "T")
    micd = load(ofile, "micd")
    its = load(ofile, "its")
    println("Data read in")

    # Convert `its` to Vector{Float64}
    its_vector = collect(its)

    # Find C from a function
    C = imm_merge_data(ps, traj, T, micd, its_vector)
    println("Data merged")

    # Check if directory exists and if not make it
    if ~isdir("Output/Imm_plots")
        mkdir("Output/Imm_plots")
    end
    # Find total number of strains
    totN = length(micd)
    #println(totN)
    # Set default plotting options
    default(dpi = 200)
    # # Plot all the populations
    # p1 = plot(yaxis = :log10,
    #           ylabel = "Population (# cells)",
    #           ylims = (1e-5, Inf),
    #           xlabel = "Time (s)")
    # for i in 1:totN
    #     # Find and eliminate zeros so that they can be plotted on a log plot
    #     inds = (C[:, i] .> 0) .& (T .<= Tmax)
    #     plot!(p1, T[inds], C[inds, i], label = "")
    # end
    # savefig(p1, "Output/Imm_plots/$(num_immigrations)events$(num_immigrants)immigrants.png")
    
    # Count number of strains at each time point
    num_strains = sum(C .> 0, dims=2)[:, 1]

    # Plot number of strains over time
    p2 = plot(T, num_strains,
              xlabel = "Time (s)",
              ylabel = "Number of Strains",
              xlims = (0, total_time/10),
              ylims = (0, maximum(num_strains) + 1),
              legend = false,
              title = "$(num_immigrations) Events, $(num_immigrants) Immigrants")
    
    savefig(p2, "Output/Imm_plots/$(num_immigrations)events$(num_immigrants)immigrants.png")
    return (nothing)
end

function multiple_freqs()
    load("Output/Imm_plots/number_of_strains.jld")
    # number of freqs 
    
    # Set maximum time to plot to
    Tmax = 5e5
    # Set default plotting options
    default(dpi = 200)
    # Plot number of strains over time
    p2 = plot(xaxis = T, y = :log10,
              xlabel = "Time (s)",
              ylabel = "Number of Strains",
              xlims = (0, total_time/10),
              ylims = (0, maximum(num_strains) + 1))
    for i in 1:totN
        # Find and eliminate zeros so that they can be plotted on a log plot
        inds = (C[:, i] .> 0) .& (T .<= Tmax)
        plot!(p2, T[inds], C[inds, i], label = "")
    end
    savefig(p2, "Output/Imm_plots/$(num_immigrations)events$(num_immigrants)immigrants.png")
    return (nothing)
end


@time frequency_plots()
@time multiple_freqs()
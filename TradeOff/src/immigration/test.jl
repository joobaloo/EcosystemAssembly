using Plots

function frequency_plots()
    # Preallocate the variables I want to extract from the input
    rps = 0
    sim_type = 0
    total_time = 6.3e7 # this is approx 2 years in seconds

    # Check that all arguments can be converted to integers
    try
        rps = parse(Int64, ARGS[1])
        sim_type = parse(Int64, ARGS[2])
    catch e
        error("Need to provide 2 integers")
    end

    if sim_type != 1
        error("Only use simulation type 1")
    end

    println("Compiled and input read in!")
    flush(stdout)
    
    # Extract other simulation parameters from the function
    Np, Nt, M, d, μrange = imm_sim_paras(sim_type)

    # Define different immigration regimes to loop through
    regimes = [
        (num_immigrations = 10, num_immigrants = 100),
        (num_immigrations = 20, num_immigrants = 200),
        (num_immigrations = 30, num_immigrants = 300)
    ]

    # Loop through each regime
    for regime in regimes
        num_immigrations = regime.num_immigrations
        num_immigrants = regime.num_immigrants

        # Read in appropriate files for multiple runs
        for r in 1:rps
            pfile = "Output/$(num_immigrations)events_$(num_immigrants)immigrants/Parameters_$r.jld"
            if !isfile(pfile)
                error("$(num_immigrants) immigrants run $r is missing a parameter file")
            end
            
            ofile = "Output/$(num_immigrations)events_$(num_immigrants)immigrants/Run$r Data.jld"
            if !isfile(ofile)
                error("$(num_immigrations) immigration events with $(num_immigrants) immigrants run $r is missing an output file")
            end

            # Read in relevant data
            ps = load(pfile, "ps")
            traj = load(ofile, "traj")
            T = load(ofile, "T")
            micd = load(ofile, "micd")
            its = load(ofile, "its")
            println("Data read in for $(num_immigrations) events and $(num_immigrants) immigrants, run $r")

            # Convert `its` to Vector{Float64}
            its_vector = collect(its)

            # Find C from a function
            C = imm_merge_data(ps, traj, T, micd, its_vector)
            println("Data merged for $(num_immigrations) events and $(num_immigrants) immigrants, run $r")

            # Count number of strains at each time point
            num_strains = sum(C .> 0, dims=2)[:, 1]

            # Set maximum time to plot to
            Tmax = 6.3e7

            # Set default plotting options
            default(dpi = 200)

            # Plot number of strains over time
            p1 = plot(T, num_strains,
                      xlabel = "Time (s)",
                      ylabel = "Number of Strains",
                      xlims = (0, Tmax),
                      ylims = (0, maximum(num_strains) + 1),
                      legend = false,
                      title = "$(num_immigrations) Events, $(num_immigrants) Immigrants, Run $r")
            
            savefig(p1, "Output/Imm_plots/num_strains_$(num_immigrations)events_$(num_immigrants)immigrants_run$r.png")
        end
    end

    return nothing
end

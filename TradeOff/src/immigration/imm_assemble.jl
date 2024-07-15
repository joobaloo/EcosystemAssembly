using TradeOff
using JLD
using Glob
include("simulation_functions.jl")
"""
    imm_assemble()
    Can't run 1 1 1 1

TBW
"""
function imm_assemble()

    # Check that sufficient arguments have been provided
    if length(ARGS) < 4
        error("Insufficient inputs provided")
    end

    # Preallocate the variables I want to extract from the input
    rps = 0
    sim_type = 0
    total_time = 6.3e7  # this is approx 2 years in seconds
    #total_time = 6.3e6
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

    if sim_type != 1 && sim_type != 5
        error("Only use simulation type 1 or 5")
    end

    # Starting run assumed to be 1
    Rs = 1

    # Check that number of repeats is greater than 0
    if rps < 1
        error("need to do at least 1 simulation")
    end
    
    # Now read in hard coded simulation parameters
    Np, Nt, M, d, μrange = imm_sim_paras(sim_type)
    #println(μrange)
    # # this μrange value ensures the thermodynamic aspect is "turned off"
    # μrange = 1.5e7 * (M / 25)

    println("Compiled and input read in!")
    flush(stdout)

    # Use formula to find number of reactions depending on number of metabolites
    if M < 4
        O = floor(Int64, M * (M - 1) / 2)
    else
        O = 4 * M - 10
    end

     # Preallocate container for filenames
     pls = fill("", Np)
     # Loop over number of required pools
     for i in 1:Np
         # Find all pools satisfying the condition
         pool_dir = joinpath(pwd(), "Pools")
         flnms = glob("ID=*N=$(Nt)M=$(M)d=$(d)u=$(μrange).jld", pool_dir)
         # Loop over valid filenames
         for j in eachindex(flnms)
             # Save first that hasn't already been used
             if flnms[j] ∉ pls
                 pls[i] = flnms[j]
             end
         end
     end
    # Save the reaction set for the first file as a point of comparison
    rs = load(pls[1], "reacs")
    # Check that all pools match this
    for i in 2:Np
        rst = load(pls[i], "reacs")
        if rst != rs
            error("pool $i uses different reaction set")
        end
    end

    # Hardcoded parameters:
    ϕR0 = 0.128 # Initial ribosome fraction is taken from ATP fits I did a while ago
    # Fairly arbitrary initial conditions
    pop = 1000.0
    conc = 1e-15
    as = 1e5
    ϕs = ϕR0
    # Starting with 10 strains for now
    Ni = 1
    
    # Time between immigration events
    mT = total_time / num_immigrations

    # Make parameter set
    ps = initialise(M, O, μrange)

    # Check that reaction set is identical to sets the pool was generated with
    if ps.reacs != rs
        error("simulation reaction set does not match pool reaction set")
    end
     
    # check directory exists
    data_dir = joinpath(
    pwd(), "Output", "$(num_immigrations)events_$(num_immigrants)immigrants")
    mkpath(data_dir)

    # Save this parameter set
    jldopen(joinpath(data_dir, "Parameters.jld"), "w") do file
        write(file, "ps", ps)
    end

    # Loading pool
    mpl = load(pls[1], "mics")
    # Now loop over the number of repeats
    for i in Rs:rps
        # Print that the new run has been started
        println("Run $i started!")
        flush(stdout)
        # Find starting time
        ti = time()
        # Then run the simulation using imm_full_simulate
        traj, T, micd, its = imm_full_simulate(ps, pop, conc, as, ϕs, mpl, Ni, mT, total_time, num_immigrations, num_immigrants)
        
        # And then print time elapsed
        tf = time()
        println("Time elapsed on run $i: $(tf - ti) s")
        # Now just save the relevant data
        jldopen(joinpath(data_dir, "Run$(i)Data.jld"), "w") do file
            # Save full set of microbe data
            write(file, "micd", micd)
            # Save extinction times
            write(file, "its", its)
            # Save time data and dynamics data
            write(file, "T", T)
            write(file, "traj", traj)
        end
        # Print to show that run has been successfully completed
        println("Run $i completed and saved!")
        flush(stdout)
    end

end

@time imm_assemble()

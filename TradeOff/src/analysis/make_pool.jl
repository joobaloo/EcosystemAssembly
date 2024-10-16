using TradeOff
export make_spool
# Function to generate a species pool
function make_spool()
    # Check that sufficient arguments have been provided
    if length(ARGS) < 3
        error("Insufficient inputs provided")
    end
    # Preallocate the variables I want to extract from the input
    #Rl = 0
    Rl = 1
    #Ru = 0
    Ru = 7
    sim_type = 0
    # Check that all arguments can be converted to integers
    try
        #Rl = parse(Int64, ARGS[1])
        #Ru = parse(Int64, ARGS[2])
        sim_type = parse(Int64, ARGS[3])
    catch e
        error("need to provide 3 integers")
    end
    # Check that simulation type is valid
    if Rl < 1
        error("lower bound on the number of reactions must be greater than 1")
    end
    # Check that simulation type is valid
    if Ru < Rl
        error("upper bound on the number of reactions can't be smaller than the lower")
    end
    # Check that a valid simulation type has been provided
    if sim_type < 1 || sim_type > 4
        error("only four simulation types defined")
    end
    # Make desired vector of reactions
    #Rs = collect(Rl:2:Ru)
    Rs = collect(Rl:Ru)
    # Extract other parameters based on simulation type chosen
    Np, Nt, M, d, μrange = sim_paras(sim_type)
    # Product to substrate ratio for equilibrium (fixing this across all simulations for now)
    mratio = 1e-2
    # Finally generate the new pool
    new_pool(Nt, M, Rs, d, μrange, mratio)
    return (nothing)
end

@time make_spool()

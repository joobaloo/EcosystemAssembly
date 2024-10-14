using TradeOff
using Random
include("../immigration/gen_pool.jl")
export make_spool
# Function to generate a species pool
function make_spool()
    #Check that sufficient arguments have been provided
    if length(ARGS) < 3
        error("Insufficient inputs provided")
    end
    #Preallocate the variables I want to extract from the input
    rl = 0
    ru = 0
    sim_type = 0

    # # Check that all arguments can be converted to integers
    try
        rl = parse(Int64, ARGS[1])
        ru = parse(Int64, ARGS[2])
        sim_type = parse(Int64, ARGS[3])
    catch e
        error("need to provide 3 integers")
    end

    # Check that simulation type is valid
    if rl < 1
        error("lower bound on the number of reactions must be greater than 1")
    end
    if ru < rl
        error("upper bound on the number of reactions can't be smaller than the lower")
    end
    if sim_type < 1 || sim_type > 4
        error("only four simulation types defined")
    end

    # Make desired vector of reactions
    Rs = collect(rl:ru)
    # Extract other parameters based on simulation type chosen
    Np, Nt, M, d, μrange = sim_paras(sim_type)
    # Product to substrate ratio for equilibrium (fixing this across all simulations for now)
    mratio = 1e-2
    # Finally generate the new pool
    μrange = 1.5e7 * (M / 25)
    imm_new_pool(Nt, M, Rs, d, μrange, mratio, rl, ru)
    println(μrange)
    return (nothing)
end

@time make_spool()

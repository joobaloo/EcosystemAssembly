using TradeOff
using JLD2
using Glob
using Plots
using Statistics

include("../immigration/simulation_functions.jl")
include("../immigration/EUE.jl")



"""
    imm_assemble()

    This function assembles the microbe community with a specific immigration rate.

    Inputs:
        This function takes 4 arguments from the command line:
            1. rps = number of repeat simulations
            2. sim_type = simulation type
            3. num_immigrations = immigration rate
            4. num_immigrants = number of immigration strains per immigration event
            5. rl = lower bound of niche size
            6. ru = upper bound of niche size

        This function also requires parameter files
    
    Outputs:
        This function produces Run#Data.jld files for each repeat of the simulation

    Note: 
        Simulation length is hardcoded to last a year (3.15e7 seconds).
"""
function imm_assemble()
    # Check that sufficient arguments have been provided
    if length(ARGS) < 4
        error("insufficient inputs provided")
    end

    # Preallocate the variables I want to extract from the input
    rps = 0
    sim_type = 0
    num_immigrations = 0
    num_immigrants = 0
    rl = 0 
    ru = 0

    # Check that all arguments can be converted to integers
    try
        rps = parse(Int64, ARGS[1])
        sim_type = parse(Int64, ARGS[2])
        num_immigrations = parse(Int64, ARGS[3])
        num_immigrants = parse(Int64, ARGS[4])
        rl = parse(Int64, ARGS[5])
        ru = parse(Int64, ARGS[6])

    catch e
        error("need to provide 6 integers")
    end
    println("ASSEMBLING COMMUNITY")
    flush(stdout)

    # Define duration of assembly
    total_time = 3.15e7 

    # Starting run assumed to be 1
    Rs = 1
    
    # Now read in hard coded simulation parameters
    Np, Nt, M, d, μrange = imm_sim_paras(sim_type)

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
         flnms = glob("ID=*Nichewidth=$(rl)-$(ru)N=$(Nt)M=$(M)d=$(d)u=$(μrange).jld", pool_dir)
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
    pwd(), "Output", "niche_size$(rl)_$(ru)", "$(num_immigrants)immigrants", "$(num_immigrations)events")
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

"""
    shan(pops::Array{Float64, 1})

    This function calculates Shannon diversity (H)
    
    Input:
        The function requires a one-dimensional array of population abundances.
    
    Output:
        This function returns the Shannon diversity index, H for a population
"""
function shan(pops::Array{Float64, 1})
    # Set initial value
    H = 0.0
    # Make vector of relative abundances
    prbs = pops ./ sum(pops)
    # Loop over relative abundances
    for i in eachindex(prbs)
        # Calculate Shannon entropy for each point
        H -= prbs[i] * log(prbs[i])
    end
    return (H)
end

function v_over_t()
    # Preallocate the variables I want to extract from the input
    rps = 0
    sim_type = 0
    num_immigrations = 0
    num_immigrants = 0
    rl = 0 
    ru = 0
    # Check that all arguments can be converted to integers
    try
        rps = parse(Int64, ARGS[1])
        sim_type = parse(Int64, ARGS[2])
        num_immigrations = parse(Int64, ARGS[3])
        num_immigrants = parse(Int64, ARGS[4])
        rl = parse(Int64, ARGS[5])
        ru = parse(Int64, ARGS[6])

    catch e
        error("need to provide 6 integers")
    end
    println("Compiled")
    flush(stdout)

    # Load in hardcoded simulation parameters
    Np, Nt, M, d, μrange = imm_sim_paras(sim_type)
  
    data_dir = joinpath(pwd(), "Output", "niche_size$(rl)_$(ru)", "$(num_immigrants)immigrants", "$(num_immigrations)events")
    # Read in parameter file
    parameter_file = joinpath(data_dir, "Parameters.jld")
    if ~isfile(parameter_file)
        error("$(num_immigrations)events_$(num_immigrants)immigrants is missing a parameter file")
    end
    # Load parameters
    ps = load(parameter_file, "ps")
    # List of pools already loaded in
    pls = []
    # Array of array to store pools
    pools = Array{Array{Microbe, 1}, 1}(undef, 1)
    # Counter for number of reactions
    NoR = 0
    # Loop over number of repeats
    for i in 1:rps
        # Load in relevant output file
        output_file = joinpath(data_dir, "Run$(i)Data.jld")
        if ~isfile(output_file)
            error("$(num_immigrations)events_$(num_immigrants)immigrants run $(i) is missing an output file")
        end
        # Load in microbe data, and immigration times
        T = load(output_file, "T")
        traj = load(output_file, "traj")
        micd = load(output_file, "micd")
        its = load(output_file, "its")

        # Convert `its` to Vector{Float64}
        its_vector = collect(its)

        # Use to construct full trajectory C
        GC.gc()
        C = imm_merge_data(ps, traj, T, micd, its_vector)
        # Preallocate vector of microbes
        ms = Array{Microbe, 1}(undef, length(micd))
        # Loop over and find each one
        for j in eachindex(micd)
            # check for case where pool hasn't already been loaded in
            if micd[j].PID ∉ pls
                # Add new pool ID in
                pls = cat(pls, micd[j].PID, dims = 1)
                # Find name of pool
                pool_file = joinpath(pwd(), "Pools",
                    "ID=$(micd[j].PID)Nichewidth=$(rl)-$(ru)N=$(Nt)M=$(ps.M)d=$(d)u=$(μrange).jld")
                # Check if this is the first pool
                if length(pls) == 1
                    # If so save the pool
                    pools[1] = load(pool_file, "mics")
                    # Find number of reactions based on this
                    NoR = maximum(pools[1] .↦ :R)
                else
                    # Otherwise just cat it on existing vector
                    pool = load(pool_file, "mics")
                    pools = cat(pools, pool, dims = 1)
                    # Find maximum number of reactions for this pool
                    NoRt = maximum(pools[1] .↦ :R)
                    # Save if higher than old number of reactions
                    NoR = max(NoR, NoRt)
                end
            end
            # Find correct pool to read from
            ind = findfirst(x -> x == micd[j].PID, pls)
            # Use this index to find and save the correct microbe
            ms[j] = (pools[ind])[micd[j].MID]
        end

        # Save total number of strains
        total_species = length(micd)

        # Preallocate containers to store number of survivors with time
        surviving_species = Array{Int64, 1}(undef, length(T))
        viable_species = Array{Int64, 1}(undef, length(T))
        total_population = zeros(length(T))
        shannon_diversity = zeros(length(T))
        no_substrates = Array{Int64, 1}(undef, length(T))
        species_per_reac_class = Array{Int64, 2}(undef, NoR, length(T))
        viable_species_per_reac_class = Array{Int64, 2}(undef, NoR, length(T))
        average_no_reac_steps = zeros(length(T))
        average_η = zeros(length(T))
        average_ω = zeros(length(T))
        total_biomass_of_viable_species = zeros(length(T))
        average_ΔG = zeros(length(T))
        average_η_per_reac_class = zeros(NoR, length(T))
        average_KS_per_reac_class = zeros(NoR, length(T))
        species_EUEs = Array{Float64, 2}(undef, length(T), total_species)
        weighted_community_EUE = Vector{Float64}(undef, 0)
        median_community_EUE =  Vector{Float64}(undef, 0)

        
        # Make vector of indices
        ϕ_i = collect((2 * total_species + ps.M + 1):(3 * total_species + ps.M))
        # Loop over all time points
        for j in 1:length(T)
            # Find indices of surviving strains
            inds = findall(x -> x > 1e-5, C[j, 1:total_species])
            # Save number of surviving strains at each time point
            surviving_species[j] = length(inds)
            # Check if at least one species survives
            if surviving_species[j] > 0
                # Save the total population at this point
                total_population[j] = sum(C[j, inds])
                # Find diversity of this point
                shannon_diversity[j] = shan(C[j, inds])
            end
            # Find indices of "viable" species
            vinds = findall(x -> x > 1e5, C[j, 1:total_species])
            # Save number of "viable" species
            viable_species[j] = length(vinds)
            # Then also number of substrates (skipping final waste product)
            no_substrates[j] = count(x -> x > 1e-12,
                C[j, (total_species + 1):(total_species + ps.M - 1)])
            # Loop over number of reactions
            for k in 1:NoR
                # Count number of strains with reaction for each case
                species_per_reac_class[k, j] = count(x -> x == k, ms[inds] .↦ :R)
                viable_species_per_reac_class[k, j] = count(x -> x == k, ms[vinds] .↦ :R)
            end
            # Set up counters for the number of species with each reaction gap
            c = zeros(M - 1)

            # loop over species
            for k in eachindex(inds)
                #empty vectors to store EUE values
                numerator_list = []
                denominator_list = []

                # loop over the reactions this strain has
                for l in 1:(ms[inds[k]].R)
                    # Find reaction number
                    Rn = ms[inds[k]].Reacs[l]
                    # Find relevant reaction
                    r = ps.reacs[Rn]

                    # println(r.Rct)
                    # println(r.Prd)
                    # println(r)
                   # println(ps)
                    # println(r.ΔG0)
                    # println(total_species)

                    # Calculate the species free-energy dissipation rate
                    D = calculate_D(
                        r.ΔG0,
                        C[j, total_species + r.Rct], 
                        C[j, total_species + r.Prd], 
                        ms[inds[k]].η[l]
                        )
                    
                    # Calculate the change in free energy of the reaction
                    ΔGT = calculate_ΔGT(
                        r.ΔG0,
                        C[j, total_species + r.Rct],
                        C[j, total_species + r.Prd],
                        )

                    # Calculate enzyme copy number
                    ECN = Eα(C[j, ϕ_i[inds[k]]], ms[inds[k]], l)
                    
                    # Calculate reaction rate
                    q = calculate_q(
                        C[j, total_species + r.Rct],
                        C[j, total_species + r.Prd],
                        ECN,
                        l,
                        ms[inds[k]],
                        r
                    )

                    numerator = (1-D/ΔGT) * q
                    denominator = q

                    push!(numerator_list, numerator)
                    push!(denominator_list, denominator)
                end
                EUE = sum(numerator_list)/sum(denominator_list)
                species_EUEs[j, inds[k]] = EUE
            end

        
            
            # Find (weighted) total eta value for viable species
            for k in eachindex(vinds)
                average_η[j] += sum(ms[vinds[k]].η .* ms[vinds[k]].ϕP) * C[j, vinds[k]]
                average_ω[j] += ms[vinds[k]].ω * C[j, vinds[k]]
                total_biomass_of_viable_species[j] += C[j, vinds[k]]
                # Loop over all reactions this strain has
                for l in 1:(ms[vinds[k]].R)
                    # Find reaction number
                    Rn = ms[vinds[k]].Reacs[l]
                    # Find relevant reaction
                    r = ps.reacs[Rn]
                    # Then calculate frac transduced
                    average_ΔG[j] += C[j, vinds[k]] * ms[vinds[k]].η[l] .*
                                     ms[vinds[k]].ϕP[l] * ΔGATP / (-r.ΔG0)
                    # Find step size
                    s_size = (r.Prd - r.Rct)
                    # weight this step size to 1 and add to total
                    average_no_reac_steps[j] += C[j, vinds[k]] * (s_size) *
                                                ms[vinds[k]].ϕP[l]
                end
            end


            # Average over biomass of viable species
            if viable_species[j] > 0
                average_η[j] /= total_biomass_of_viable_species[j]
                average_ω[j] /= total_biomass_of_viable_species[j]
                average_no_reac_steps[j] /= total_biomass_of_viable_species[j]
                average_ΔG[j] /= total_biomass_of_viable_species[j]
            end
            # Break down eta and omega value by R
            for k in eachindex(vinds)
                # Find relevant reaction number
                l = ms[vinds[k]].R
                # Add contribution to relevant total
                average_η_per_reac_class[l, j] += sum(ms[vinds[k]].η .* ms[vinds[k]].ϕP)
                average_KS_per_reac_class[l, j] += sum(ms[vinds[k]].KS .* ms[vinds[k]].ϕP)
            end
            # Now weight by number of strains with each type of reaction
            for k in 1:NoR
                if viable_species_per_reac_class[k, j] > 0
                    average_η_per_reac_class[k, j] /= viable_species_per_reac_class[k, j]
                    average_KS_per_reac_class[k, j] /= viable_species_per_reac_class[k, j]
                end
            end
            
        end

       

        # Preallocate final ϕR values
        final_ϕR = zeros(surviving_species[end])
        # Find indices of all surviving species
        inds = findall(x -> x > 1e-5, C[end, 1:total_species])
        # Loop over number of survivors
        for j in 1:surviving_species[end]
            # Store corresponding final ϕR value
            final_ϕR[j] = C[end, ϕ_i[inds[j]]]
        end

        # Find mean Community EUE
        for j in 1:length(T)
            community_EUE = []
            # Find indices of surviving strains
            inds = findall(x -> x > 1e-5, C[j, 1:total_species])
            for k in inds
                if !isnan(C[j,k])
                    species_EUE = (C[j, k] * species_EUEs[j, k])/sum(filter(!isnan, (C[j, inds])))
                    push!(community_EUE, species_EUE)
                end
            end
            if any(!isnan, community_EUE)
                push!(weighted_community_EUE, sum(filter(!isnan, community_EUE)))
            else
                push!(weighted_community_EUE, 0.0)
            end
        end

        # find median community_EUE
        for j in 1:length(T)

            # Find indices of surviving strains
            inds = findall(x -> x > 1e-5, C[j, 1:total_species])

            if all(isnan, species_EUEs[j, inds])
                push!(median_community_EUE, 0.0)
            else 
                push!(median_community_EUE, median(filter(!isnan,species_EUEs[j, inds])))
            end
        end

        # Now just save the relevant data
        jldopen(joinpath(data_dir, "AvRun$(i)Data.jld"), "w") do file
            # Save full time course
            write(file, "T", T)
            # Save reaction data
            write(file, "species_per_reac_class", species_per_reac_class)
            write(file, "viable_species_per_reac_class", viable_species_per_reac_class)
            write(file, "average_η_per_reac_class", average_η_per_reac_class)
            write(file, "average_KS_per_reac_class", average_KS_per_reac_class)
            # Save the other quantities
            write(file, "surviving_species", surviving_species)
            write(file, "viable_species", viable_species)
            write(file, "total_population", total_population)
            write(file, "total_biomass_of_viable_species", total_biomass_of_viable_species)
            write(file, "shannon_diversity", shannon_diversity)
            write(file, "no_substrates", no_substrates)
            write(file, "average_η", average_η)
            write(file, "average_no_reac_steps", average_no_reac_steps)
            write(file, "average_ω", average_ω)
            write(file, "average_ΔG", average_ΔG)
            write(file, "final_ϕR", final_ϕR)
            # Save EUE data
            write(file, "species_EUEs", species_EUEs)
            write(file, "community_EUE", weighted_community_EUE)
            write(file, "median_community_EUE", median_community_EUE)
            # Finally save final time to help with benchmarking
            write(file, "final_time_point", T[end])
        end
        println("Run $i analysed")
        flush(stdout)
    end
    return (nothing)
end

# Function to make a dictionary to store all desired data based on list of names and
# dimensions. This also involves preallocating data where relevant
function make_data_dictonary_from_list(
        variables_of_interest::Array{
            Tuple{String, Int64,
                Bool, Bool}, 1
        },
        repeats::Int64,
        no_reactions::Int64,
        no_time_points::Int64)
    # Make initially empty dictionary
    data_dict = Dict()
    # Then use a for loop to populate it with names + preallocated memory
    for variable in variables_of_interest
        if variable[2] == 2
            data_dict[variable[1]] = Dict(
                "combined_data" => zeros(repeats, no_time_points),
                "run_data" => Float64[], "dims" => 2,
                "viable" => variable[3],
                "divide_by_R" => variable[4])
        elseif variable[2] == 3
            data_dict[variable[1]] = Dict(
                "combined_data" => zeros(repeats, no_reactions,
                    no_time_points),
                "run_data" => Float64[], "dims" => 3,
                "viable" => variable[3],
                "divide_by_R" => variable[4])
        end
    end

    return data_dict
end

# This function loads in the relevant data and adds it into the data dictionary
function load_trajectory_vars_to_dict!(vfile::String, data_dict::Dict)
    # Loop over every name in the data dictionary
    for variable in keys(data_dict)
        data_dict[variable]["run_data"] = load(vfile, variable)
    end
    return (data_dict)
end

# This functions uses interpolation to take time slices for each variable. This is done to
# ensure that runs can be sensibly compared
function add_to_combined_data!(data_dict::Dict, T::Array{Float64, 1},
        times::Array{Float64, 1}, Tind::Int64, cnt::Int64, i::Int64)
    # Skip averaging if previous point is missing
    if Tind > 1 && data_dict["viable_species"]["run_data"][Tind - 1] != 0
        # Calculate relevant time gaps
        Tg = (T[Tind] - T[Tind - 1])
        T1x = times[cnt] - T[Tind - 1]
        T2x = T[Tind] - times[cnt]
        # Loop over every name in the data dictionary
        for variable in keys(data_dict)
            # Check variable dimensionality
            if data_dict[variable]["dims"] == 2
                # Interpolate and add into combined data
                point = interpolate_time(data_dict[variable]["run_data"], Tg, T1x, T2x,
                    Tind)
                data_dict[variable]["combined_data"][i, cnt] = point
            elseif data_dict[variable]["dims"] == 3
                point = interpolate_time(data_dict[variable]["run_data"], Tg, T1x, T2x,
                    Tind)
                data_dict[variable]["combined_data"][i, :, cnt] = point
            end
        end
    else
        # Loop over every name in the data dictionary
        for variable in keys(data_dict)
            # Check variable dimensionality
            if data_dict[variable]["dims"] == 2
                # In the one case just add the value at time = 0
                point = data_dict[variable]["run_data"][Tind]
                data_dict[variable]["combined_data"][i, cnt] = point
            elseif data_dict[variable]["dims"] == 3
                point = data_dict[variable]["run_data"][:, Tind]
                data_dict[variable]["combined_data"][i, :, cnt] = point
            end
        end
    end
    return (data_dict)
end

# This function calculates the means for all variables along the trajectory
function calculate_trajectory_means!(data_dict::Dict, no_reactions::Int64,
        no_simulations::Vector{Float64},
        no_viable_simulations::Vector{Float64},
        no_simulations_with_R::Array{Float64, 2})
    # Loop over every name in the data dictionary
    for variable in keys(data_dict)
        # Sum to find totals
        total = dropdims(sum(data_dict[variable]["combined_data"], dims = 1), dims = 1)
        # Check variable dimensionality (3D case requires preallocation)
        if data_dict[variable]["dims"] == 2
            # Now calculate means
            if data_dict[variable]["viable"]
                data_dict[variable]["means"] = total ./ no_viable_simulations
            else
                data_dict[variable]["means"] = total ./ no_simulations
            end
        elseif data_dict[variable]["dims"] == 3
            means = zeros(no_reactions, size(total, 2))
            # Calculate means (first checking what should be divided by)
            if data_dict[variable]["divide_by_R"]
                for i in 1:no_reactions
                    means[i, :] = total[i, :] ./ no_simulations_with_R[i, :]
                end
            elseif data_dict[variable]["viable"]
                for i in 1:no_reactions
                    means[i, :] = total[i, :] ./ no_viable_simulations
                end
            else
                for i in 1:no_reactions
                    means[i, :] = total[i, :] ./ no_simulations
                end
            end
            # Add means to data dictionary
            data_dict[variable]["means"] = means
        end
    end
    return (data_dict)
end

# Written own function for standard deviation as it was tricky to get the standard function
# to handle no species simulations properly
function find_stand_dev(values::Vector{Float64}, mean::Float64, no_points::Float64)
    return sqrt(sum((values .- mean) .^ 2) / (no_points - 1))
end

# This function calculates the standard deviations for all variables along the trajectory
function calculate_trajectory_standard_devs!(data_dict::Dict, times::Vector{Float64},
        final_time_points::Vector{Float64},
        no_simulations::Vector{Float64},
        no_viable_simulations::Vector{Float64},
        no_simulations_with_R::Matrix{Float64},
        no_reactions::Int64)
    # Preallocate containers for the standard deviations
    for variable in keys(data_dict)
        data_dict[variable]["sds"] = zeros(size(data_dict[variable]["means"]))
    end
    # Loop over times
    for i in eachindex(times)
        # Find indices of still progressing trajectories
        inds = (final_time_points .>= times[i])
        # Find indices of still progressing trajectories with one or more viable strains
        vinds = (final_time_points .>= times[i]) .&
                (data_dict["viable_species"]["combined_data"][:, i] .> 0.0)
        # Find indices for where viable species in each reaction class exist
        rinds = Vector{Vector}(undef, no_reactions)
        for j in 1:no_reactions
            rinds[j] = (final_time_points .>= times[i]) .&
                       (data_dict["viable_species_per_reac_class"]["combined_data"][:, j,
                i] .>
                        0.0)
        end
        # Loop over all variables
        for variable in keys(data_dict)
            # Calculate standard deviations (procedure changes based on what needs to be
            # averaged over)
            if data_dict[variable]["dims"] == 2 && data_dict[variable]["viable"]
                # These should be calculated just for viable strains
                if no_viable_simulations[i] > 1
                    sd_value = find_stand_dev(
                        data_dict[variable]["combined_data"][vinds,
                            i],
                        data_dict[variable]["means"][i],
                        no_viable_simulations[i])
                    data_dict[variable]["sds"][i] = sd_value
                else
                    data_dict[variable]["sds"][i] = NaN
                end
            elseif data_dict[variable]["dims"] == 2
                sd_value = find_stand_dev(data_dict[variable]["combined_data"][inds, i],
                    data_dict[variable]["means"][i],
                    no_simulations[i])
                data_dict[variable]["sds"][i] = sd_value
            elseif data_dict[variable]["dims"] == 3 && data_dict[variable]["divide_by_R"]
                # Calculate standard deviations for reactions
                for j in 1:no_reactions
                    # Use only these in the reaction calculation
                    if no_simulations_with_R[j, i] > 1
                        sd_value = find_stand_dev(
                            data_dict[variable]["combined_data"][rinds[j],
                                j,
                                i],
                            data_dict[variable]["means"][j, i],
                            no_simulations_with_R[j, i])
                        data_dict[variable]["sds"][j, i] = sd_value
                    else
                        data_dict[variable]["sds"][j, i] = NaN
                    end
                end
            elseif data_dict[variable]["dims"] == 3 && data_dict[variable]["viable"]
                # These should be calculated just for viable strains
                if no_viable_simulations[i] > 1
                    for j in 1:no_reactions
                        sd_value = find_stand_dev(
                            data_dict[variable]["combined_data"][vinds,
                                j,
                                i],
                            data_dict[variable]["means"][j, i],
                            no_viable_simulations[i])
                        data_dict[variable]["sds"][j, i] = sd_value
                    end
                else
                    data_dict[variable]["sds"][:, i] .= NaN
                end
            elseif data_dict[variable]["dims"] == 3
                for j in 1:no_reactions
                    sd_value = find_stand_dev(
                        data_dict[variable]["combined_data"][inds, j,
                            i],
                        data_dict[variable]["means"][j, i],
                        no_simulations[i])
                    data_dict[variable]["sds"][j, i] = sd_value
                end
            end
        end
    end
    return (data_dict)
end

# Function to read in variables with time and calculate stats over time
"""
    calculate_trajectory_stats()

TBW
"""
function calculate_trajectory_stats()
    # Check that sufficient arguments have been provided
    if length(ARGS) < 4
        error("insufficient inputs provided")
    end

    # Preallocate the variables I want to extract from the input
    repeats = 0
    sim_type = 0
    num_immigrations = 0
    num_immigrants = 0
    rl = 0 
    ru = 0
    # Check that all arguments can be converted to integers
    try
        repeats = parse(Int64, ARGS[1])
        sim_type = parse(Int64, ARGS[2])
        num_immigrations = parse(Int64, ARGS[3])
        num_immigrants = parse(Int64, ARGS[4])
        rl = parse(Int64, ARGS[5])
        ru = parse(Int64, ARGS[6])

    catch e
        error("need to provide 6 integers")
    end

     println("Compiled")
     flush(stdout)
    
    # Load in hardcoded simulation parameters
    Np, Nt, M, d, μrange = sim_paras(sim_type)

    # Number of steps to calculate stats for
    no_steps = 2500

    # Container to store final times
    final_time_points = zeros(repeats)

    # Define here so that it is available outside the for loop
    no_reactions = 0

    # Define data directory
    data_dir = joinpath(pwd(), "Output", "niche_size$(rl)_$(ru)", "$(num_immigrants)immigrants", "$(num_immigrations)events")
    # Loop over number of repeats
    for i in 1:repeats
        # Load in relevant output file
        averages_file = joinpath(data_dir, "AvRun$(i)Data.jld")
        if ~isfile(averages_file)
            error("$(num_immigrations)events_$(num_immigrants)immigrants run $(i) is missing a variables file")
        end
        # Just want to save final times for now
        final_time_points[i] = load(averages_file, "final_time_point")
        # Save number of reactions from the first run
        if i == 1
            viable_species_per_reac_class = load(
                averages_file, "viable_species_per_reac_class")
            no_reactions = size(viable_species_per_reac_class, 1)
        end
    end
    # Use maximum final time to set time value
    times = collect(range(0.0, maximum(final_time_points), length = no_steps))
    # Preallocate relevant containers
    no_simulations = zeros(length(times))
    no_viable_simulations = zeros(length(times))
    no_simulations_with_R = zeros(no_reactions, length(times))
    # Define the other variables of interest + their required dimensionality
    # + whether to divide by only simulations with viable species
    # + whether to divide by number of species in reaction number class
    variables_of_interest = [
        ("surviving_species", 2, false, false),
        ("viable_species", 2, false, false),
        ("total_population", 2, false, false),
        ("shannon_diversity", 2, false, false),
        ("no_substrates", 2, false, false),
        ("species_per_reac_class", 3, false, false),
        ("viable_species_per_reac_class", 3, true, false),
        ("average_no_reac_steps", 2, true, false),
        ("average_η", 2, true, false),
        ("average_ω", 2, true, false),
        ("total_biomass_of_viable_species", 2, true, false),
        ("average_ΔG", 2, true, false),
        ("average_η_per_reac_class", 3, true, true),
        ("average_KS_per_reac_class", 3, true, true),
        ("community_EUE", 2, false, false),
        ("median_community_EUE", 2, false, false)
    ]
    # Convert this list into a dictionary of preallocated arrays
    data_dict = make_data_dictonary_from_list(variables_of_interest, repeats,
        no_reactions, length(times))
    # Final ϕR values are a special case
    all_final_ϕRs = Float64[]
    # Loop over number of trajectories (to minimise the number of reads in)
    for i in 1:repeats
        # Load in relevant output file
        averages_file = joinpath(data_dir, "AvRun$(i)Data.jld")
        if ~isfile(averages_file)
            error("$(num_immigrations)events_$(num_immigrants)immigrants run $(rN) is missing a variables file")
        end
        # First find and save final ϕR value for run
        final_ϕR = load(averages_file, "final_ϕR")
        all_final_ϕRs = cat(all_final_ϕRs, final_ϕR, dims = 1)
        # Then load the time interval data
        T = load(averages_file, "T")
        # Load all other variables into the dictionary
        data_dict = load_trajectory_vars_to_dict!(averages_file, data_dict)
        # Bool to indicate end of the run
        run_end = false
        cnt = 0
        # Loop until end of run reached
        while run_end == false
            # increment counter
            cnt += 1
            # Still a trajectory here so increment that count by one
            no_simulations[cnt] += 1
            # Find index of first point greater than or equal to time
            Tind = findfirst(x -> x >= times[cnt], T)
            # Only increment the viable counter if there are viable strains at time point
            if data_dict["viable_species"]["run_data"][Tind] != 0
                no_viable_simulations[cnt] += 1
            end
            # Loop over reactions to find number of viable reactions
            for j in 1:no_reactions
                # Check if counter should be incremented
                if data_dict["viable_species_per_reac_class"]["run_data"][j, Tind] .> 0.0
                    no_simulations_with_R[j, cnt] += 1
                end
            end
            # Use function to add data into combined data container
            data_dict = add_to_combined_data!(data_dict, T, times, Tind, cnt, i)
            # Finally check if next time point is higher than final time for this trajectory
            if cnt >= length(times) || times[cnt + 1] > final_time_points[i]
                run_end = true
            end
        end
        println("Analysed trajectory $(i)")
        flush(stdout)
    end
    # Use function to calculate means
    data_dict = calculate_trajectory_means!(data_dict, no_reactions, no_simulations,
        no_viable_simulations, no_simulations_with_R)
    println("Means found")
    # Then use function to calculate standard deviations
    data_dict = calculate_trajectory_standard_devs!(data_dict, times, final_time_points,
        no_simulations, no_viable_simulations,
        no_simulations_with_R, no_reactions)
    println("Standard deviations found")
    # Now want to save means and standard deviations
    jldopen(joinpath(data_dir, "RunStats$(num_immigrations)events_$(num_immigrants)immigrants.jld"), "w") do file
        # Save times
        write(file, "times", times)
        # Save number of continuing trajectories
        write(file, "no_simulations", no_simulations)
        write(file, "no_viable_simulations", no_viable_simulations)
        write(file, "no_simulations_with_R", no_simulations_with_R)
        # Save means and standard deviations
        for variable in keys(data_dict)
            write(file, "mean_$(variable)", data_dict[variable]["means"])
            write(file, "sd_$(variable)", data_dict[variable]["sds"])
        end
        # Finally write all of the final ϕR values out
        write(file, "all_final_ϕRs", all_final_ϕRs)
    end
    println("All data saved")
    return (nothing)
end


@time begin
    imm_assemble()
    v_over_t()
    calculate_trajectory_stats()
end 

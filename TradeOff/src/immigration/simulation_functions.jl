using DifferentialEquations
include("../simulate/sim.jl")

export initialise, choose_reactions, imm_merge_data

function imm_merge_data(ps::TOParameters,
    traj::Array{Array{Float64, 2}, 1},
    T::Array{Float64, 1},
    micd::Array{MicData, 1},
    its::Array{Float64, 1})
    # Find total number of microbes
    totN = length(micd)
    # Find total number of immigration attempts
    ims = length(its)
    # Preallocate array to store all the trajectory data
    C = Array{Float64, 2}(undef, length(T), 3 * totN + ps.M)
    # Previous immigration time is initially zero
    tp = 0.0
    # Index is 1 to begin with
    ind_tp = 1
    # Then loop over every immigration attempt
    for i in 1:ims
        # Extract relevant trajectory
        tt = traj[i]
        # Find new immigration time
        tn = its[i]
        # Find index of this time in vector T
        ind_tn = ind_tp + size(tt, 1) - 2
        # Find strains that exist in this window
        inds = ((micd .↦ :ImT) .<= tp) .& (((micd .↦ :ExT) .>= tn) .| isnan.(micd .↦ :ExT))
        # From this calculate number of strains
        Ns = sum(inds)
        # setup counter
        cnt = 1
        # Firstly save the concentrations
        C[ind_tp:ind_tn, (totN + 1):(totN + ps.M)] = tt[1:(end - 1), (Ns + 1):(Ns + ps.M)]
        # loop over total number of strains
        for j in 1:totN
            # Find strains that exist within this window
            if inds[j] == true
                C[ind_tp:ind_tn, j] = tt[1:(end - 1), cnt]
                C[ind_tp:ind_tn, totN + ps.M + j] = tt[1:(end - 1), Ns + ps.M + cnt]
                C[ind_tp:ind_tn, 2 * totN + ps.M + j] = tt[1:(end - 1), 2 * Ns + ps.M + cnt]
                # Then increment counter
                cnt += 1
            else
                # Strains that aren't present have their variables set as NaN
                C[ind_tp:ind_tn, j] .= NaN
                C[ind_tp:ind_tn, totN + ps.M + j] .= NaN
                C[ind_tp:ind_tn, 2 * totN + ps.M + j] .= NaN
            end
        end
        # Finally update previous time
        tp = tn
        # Update previous index to be one higher than final index last time
        ind_tp = ind_tn + 1
    end
    # find strains that go extinct at the very end or survive past it
    inds = (((micd .↦ :ExT) .== T[end]) .| isnan.(micd .↦ :ExT))
    # Save number of survivors
    Ns = sum(inds)
    # Extract final relevant trajectory
    tt = traj[ims + 1]
    # setup counter
    cnt = 1
    # Firstly save the concentrations
    C[ind_tp:end, (totN + 1):(totN + ps.M)] = tt[1:end, (Ns + 1):(Ns + ps.M)]
    # loop over total number of strains
    for j in 1:totN
        # Find strains that exist within this window
        if inds[j] == true
            C[ind_tp:end, j] = tt[1:end, cnt]
            C[ind_tp:end, totN + ps.M + j] = tt[1:end, Ns + ps.M + cnt]
            C[ind_tp:end, 2 * totN + ps.M + j] = tt[1:end, 2 * Ns + ps.M + cnt]
            # Then increment counter
            cnt += 1
        else
            # Strains that aren't present have their variables set as NaN
            C[ind_tp:end, j] .= NaN
            C[ind_tp:end, totN + ps.M + j] .= NaN
            C[ind_tp:end, 2 * totN + ps.M + j] .= NaN
        end
    end
    return (C)
end

# function to randomly choose the reactions that a specific microbe can make use of
function choose_reactions(O::Int64, Rs::Array{Int64, 1})
    @assert length(Rs)>0 "vector of possible numbers of reactions can't be empty"
    # select (uniformly) R value from supplied vector
    R = rand(Rs)
    # Preallocate vector of reaction identities
    Reacs = zeros(Int64, R)
    # Choose random first value
    Reacs[1] = rand(1:O)
    # Then fill out later values
    for i in 2:R
        good = false
        while good == false
            r = rand(1:O)
            # Check to avoid repeated values
            if r ∉ Reacs[1:(i - 1)]
                Reacs[i] = r
                good = true
            end
        end
    end
    return (R, Reacs)
end

# function to generate parameter set for the fixed parameters
function initialise(M::Int64, O::Int64)
    # Assume that temperature T is constant at 20°C
    T = 293.15
    # From Posfai et al (2017) dilution rate 0.21 per hour
    δ = 6.0e-5 * ones(M) # Metabolite dilution rate
    # Human blood glucose is approximately 5.5 m mol per litre (wikipedia)
    # Sensible order of magnitude to aim for set κ/δ = 5.5e-3
    κ = zeros(M)
    # All but resource 1 is not supplied
    κ[1] = 3.3e-7 # Metabolite supply rate
    # Chosen so that 100 steps yields slightly more free energy than respiring glucose
    #μrange = 5e6 * (M / 25)
    μrange = 1.5e7 * (M / 25)
    # Generate fixed set of reactions
    RP, ΔG = fix_reactions(O, M, μrange, T)
    # Preallocate vector of reactions
    reacs = Array{Reaction, 1}(undef, O)
    for i in 1:O
        reacs[i] = make_Reaction(i, RP[i, 1], RP[i, 2], ΔG[i])
    end
    # Now make the parameter set
    ps = make_TOParameters(M, O, T, κ, δ, reacs)
    return (ps)
end

# function to generate parameter set for the fixed parameters
function initialise(M::Int64, O::Int64, μrange::Float64)
    # Assume that temperature T is constant at 20°C
    T = 293.15
    # From Posfai et al (2017) dilution rate 0.21 per hour
    δ = 6.0e-5 * ones(M) # Metabolite dilution rate
    # Human blood glucose is approximately 5.5 m mol per litre (wikipedia)
    # Sensible order of magnitude to aim for set κ/δ = 5.5e-3
    κ = zeros(M)
    # All but resource 1 is not supplied
    κ[1] = 3.3e-7 # Metabolite supply rate
    # Generate fixed set of reactions
    RP, ΔG = fix_reactions(O, M, μrange, T)
    # Preallocate vector of reactions
    reacs = Array{Reaction, 1}(undef, O)
    for i in 1:O
        reacs[i] = make_Reaction(i, RP[i, 1], RP[i, 2], ΔG[i])
    end
    # Now make the parameter set
    ps = make_TOParameters(M, O, T, κ, δ, reacs)
    return (ps)
end

"""
    imm_sim_paras(sim_type::Int64)

TBW
"""
function imm_sim_paras(sim_type::Int64)
    # Set the hardcoded variables here
    Np = 1
    Nt = 5000
    M = 25
    μrange = 1.5e7 * (M / 25)
    # Choose biomass loss rate and
    if sim_type == 1 || sim_type == 5
        d = 6e-5
        #μrange = 5e6 * (M / 25)
    elseif sim_type == 2
        d = 6e-5
        #μrange = 1.5e6 * (M / 25)
    elseif sim_type == 3
        d = 1e-4
        #μrange = 5e6 * (M / 25)
    elseif sim_type == 4
        d = 1e-4
        #μrange = 1.5e6 * (M / 25)
    end
    return (Np, Nt, M, d, μrange)
end

"""
    imm_full_simulate(ps::TOParameters, 
                  pop::Float64,
                  conc::Float64,
                  as::Float64,
                  ϕs::Float64,
                  mpl::Array{Microbe, 1},
                  Ni::Int64,
                  mT::Float64,
                  total_time::Float64,      
                  num_immigrations::Int64,  
                  num_immigrants::Int64)    

Simulates the dynamics of a microbial community with immigration events over a specified time period.

Parameters:
- `ps::TOParameters`: Parameters for the simulation dynamics.
- `pop::Float64`: Initial population size for each microbe.
- `conc::Float64`: Initial concentration of resources.
- `as::Float64`: Internal energy concentrations.
- `ϕs::Float64`: Ribosome fractions.
- `mpl::Array{Microbe, 1}`: Pool of microbes.
- `Ni::Int64`: Initial number of microbe strains in the community.
- `mT::Float64`: Average time between immigration events.
- `total_time::Float64`: Total time duration of the simulation.
- `num_immigrations::Int64`: Number of immigration events during the simulation.
- `num_immigrants::Int64`: Number of immigrants per immigration event.

Returns a tuple `(traj, T, micd, its)`:
- `traj`: Array of trajectories of microbial population dynamics.
- `T`: Time points corresponding to each trajectory.
- `micd`: Array of `MicData` capturing microbe-specific data over time.
- `its`: Array of immigration times.

The function simulates the growth and dynamics of a microbial community with periodic immigration events, evenly spread across `total_time`. Each immigration event introduces `num_immigrants` new microbe types into the community.

"""
function imm_full_simulate(ps::TOParameters, 
    pop::Float64,
    conc::Float64,
    as::Float64,
    ϕs::Float64,
    mpl::Array{Microbe, 1},
    Ni::Int64,
    mT::Float64,
    total_time::Float64, # Total simulation time
    num_immigrations::Int64, # Number of immigration events
    num_immigrants::Int64) # Number of immigrants per event

    # Make container to store microbial data
    micd = Array{MicData}(undef, Ni)

    # Make container to store trajectory data
    traj = Array{Array{Float64, 2}}(undef, num_immigrations + 1)

    # Set a value for the maximum number of strains that can be simulated
    max_N = 300

    # Preallocate memory
    rate = Array{Float64, 2}(undef, max_N, ps.O)

    # Now substitute preallocated memory in
    dyns!(dx, x, ms, t) = full_dynamics!(dx, x, ms, ps, rate, t)

    # Preallocate initial vector of microbes
    ms = Array{Microbe, 1}(undef, Ni)

    # Bool to store if full vector of microbes has been found
    full = false
    i = 1
    # Randomly choose them from the pool
    while full == false
        r = rand(1:length(mpl))

        # Check that strain hasn't already been added
        if mpl[r] ∉ ms[1:(i - 1)]
            # And if not add it
            ms[i] = mpl[r]
            # Increment counter
            i += 1
        end

        # End while loop when all the microbes have been filled
        if i == length(ms) + 1
            full = true
        end
    end

    # Store generated microbes as MicData
    for i in 1:Ni
        micd[i] = make_MicData(ms[i].ID, ms[i].PID, 0.0, NaN)
    end

    # Initial number of surviving strains is equal to 1
    Ns = Ni

    # Make initial values
    pops = pop * ones(length(ms))
    concs = conc * ones(ps.M)
    ass = as * ones(length(ms))
    ϕss = ϕs * ones(length(ms))
    x0 = [pops; concs; ass; ϕss]

    # Preallocate immigration times
    its = range(0, stop=total_time, length=num_immigrations+1)[2:end]  # Evenly spread immigration times
    
    # Define initial step
    tspan = (0, its[1])

    # Then setup and solve the initial problem
    prob = ODEProblem(dyns!, x0, tspan, ms)
    sol = DifferentialEquations.solve(prob)

    # Make containers to store dynamics
    T = sol.t
    C = copy(sol')

    # Save this C for output
    traj[1] = C[:, :]

    # Find indices of surviving strains
    inds = sol'[end, 1:Ni] .> 1e-5

    # Make vector to store indices to delete
    dls = []

    # Find any extinctions
    for i in 1:Ni
        if inds[i] == false
            # Mark extinction time in the microbe data
            micd[i] = make_MicData(micd[i].MID, micd[i].PID, micd[i].ImT, its[1])
            # Mark species for deletion
            dls = cat(dls, i, dims = 1)
            # Set extinct species values as NaN in the output data
            C[end, i] = NaN
            C[end, ps.M + Ni + i] = NaN
            C[end, ps.M + 2 * Ni + i] = NaN
            # Reduce number of surviving strains counter by 1
            Ns -= 1
        end
    end

    # If no immigration events are considered then ensure that changes to C are retained
    if num_immigrations == 0
        traj[1] = C[:, :]
    end

    # Delete extinct species
    ms = deleteat!(ms, dls)

    # Now loop over for every immigration attempt
    for i in 1:num_immigrations
        println("Immigration attempt $i initiated")
            flush(stdout)
        # if i % 50 == 0
        #     println("Immigration attempt $i initiated")
        #     flush(stdout)
        # end

        # Use fixed number of immigrants
        nI = num_immigrants

        # Make new vector to store microbes
        mst = Array{Microbe, 1}(undef, nI)

        # Make container to store microbial data
        micdt = Array{MicData}(undef, nI)

        # Set up while loop to check microbes
        j = 1
        full = false

        # Randomly choose them from the pool
        while full == false
            r = rand(1:length(mpl))

            # Check that strain hasn't already been added
            if mpl[r] ∉ ms && mpl[r] ∉ mst[1:(j - 1)]

            # And if not add it
            mst[j] = mpl[r]
            
            # Increment counter
            j += 1
        end
        # End while loop when all the microbes have been filled
        if j == length(mst) + 1
            full = true
        end
    end

    # Add microbes to the existing vector
    ms = cat(ms, mst, dims = 1)

    # Find new MicData values
    for j in 1:nI
        micdt[j] = make_MicData(mst[j].ID, mst[j].PID, its[i], NaN)
    end

    # Add this new data to the old
    micd = cat(micd, micdt, dims = 1)

    # Find time to next immigration, if not at the last step
    if i != num_immigrations
        # Use previous immigration time to define the time step
        tspan = (its[i], its[i + 1])
    else
        # At last step just integrate for five times the average time, so that dynamics settle
        tf = total_time + 5 * mT
        # Use previous immigration time to define the time span
        tspan = (its[i], tf)
    end

    # Find all indices of still relevant initial conditions in C
    in_cons = findall(!isnan, C[end, :])

    # Then find initial conditions directly from C
    pops_old = C[end, in_cons[1:Ns]]
    concs = C[end, in_cons[(Ns + 1):(Ns + ps.M)]]
    as_old = C[end, in_cons[(Ns + ps.M + 1):(2 * Ns + ps.M)]]
    ϕs_old = C[end, in_cons[(2 * Ns + ps.M + 1):(3 * Ns + ps.M)]]

    # Find indices of concentrations below threshold and overwrite
    nCids = findall(x -> x < 1e-15, concs)

    # Then set all these to zero
    concs[nCids] .= 1e-15

    # Make new vectors incorporating old and new microbes
    pops = cat(pops_old, pop * ones(length(mst)), dims = 1)
    ass = cat(as_old, as * ones(length(mst)), dims = 1)
    ϕss = cat(ϕs_old, ϕs * ones(length(mst)), dims = 1)

    # Collect all of this together in a vector of initial conditions
    x0 = [pops; concs; ass; ϕss]

    # Now setup and solve the problem with the new strains
    prob = ODEProblem(dyns!, x0, tspan, ms)
    sol = DifferentialEquations.solve(prob)

    # Update the number of survivors, as new strains have been added
    Ns += nI

    # Store new dynamics in a temporary form
    Tt = sol.t
    C = copy(sol')

    # Save new dynamics for output
    traj[i + 1] = C

    # Add to full vector of times
    T = cat(T, Tt[2:end], dims = 1)

    # Now find indices of recently extinct strains
    inds = (C[end, 1:Ns] .<= 1e-5)

    # Make vector to store indices to delete
    dls = []

    # Find indices of all strains that still survive in micd
    svs = findall(isnan, (micd .↦ :ExT))

    # Find any extinctions
    for j in 1:Ns
        if inds[j] == true
            # Find index of the newly extinct strain
            ex = svs[j]
            # The essential problem is converting j into the true index
            # Mark extinction time in the microbe data
            micd[ex] = make_MicData(micd[ex].MID, micd[ex].PID, micd[ex].ImT, tspan[2])
            # Find indices of species with ID's matching the one being made extinct
            sinds = findall(x -> x == micd[ex].MID, ms .↦ :ID)
            # Check if there's multiple
            if length(sinds) > 1
                # If there is compare pool ids
                pind = findfirst(x -> x == micd[ex].PID, ms[sinds] .↦ :PID)
                sind = sinds[pind]
            else
                sind = sinds[1]
            end
            # Mark species for deletion
            dls = cat(dls, sind, dims = 1)
            # Set extinct species values as NaN in the output data
            C[end, j] = NaN
            C[end, ps.M + Ns + j] = NaN
            C[end, ps.M + 2 * Ns + j] = NaN
        end
    end
    
    # Reduce number of surviving strains by the number of extinctions
    Ns -= sum(inds)
    # Delete extinct species
    ms = deleteat!(ms, dls)
    end
    return (traj, T, micd, its)
end

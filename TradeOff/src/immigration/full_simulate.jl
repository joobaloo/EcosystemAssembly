#Simulation code to run one instance of the simulation with a user defined starting condition
# ps is parameter set, Tmax is the time to integrate to, pop, conc, as and ϕs
# are the initial conditions, mpl is a pool of microbes, mT is mean immigration time,
# ims is the number of immigrations, λIm controls rate of additional immigrants

function full_simulate(ps::TOParameters, # parameter set
                       pop::Float64, # initial population size
                       conc::Float64, # initial concentration
                       as::Float64,  #??
                       ϕs::Float64, #??
                       mpl::Array{Microbe, 1}, # pool of microbes
                       Ni::Int64, # initial number of microbes
                       mT::Float64, # mean time between immigration events
                       ims::Int64, # number of immigration events
                       λIm::Float64,
                       sim_length::Float64) #number of invading species same as nI kinda
    # Preallocate immigration times
    its = 
    # Make container to store microbial data
    micd = Array{MicData}(undef, Ni)
    # Make container to store trajectory data
    traj = Array{Array{Float64, 2}}(undef, ims + 1)
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
    # Make distribution to sample random immigration times from
    td = Exponential(mT)
    # Make distribution to sample random number of invading species from
    sd = Poisson(λIm)
    # Check if this is a no-immigration simulation
    if ims == 0
        # In this case integrate for five times the average time, so that dynamics settle
        ti = 25 * mT
    else
        # Otherwise choose a random time for the initial step
        ti = rand(td)
        # Save this as the first immigration time
        its[1] = ti
    end
    # Define initial step
    tspan = (0, ti)
    # Then setup and solve the initial problem
    prob = ODEProblem(dyns!, x0, tspan, ms)
    sol = DifferentialEquations.solve(prob)
    # Make containers to store dynamics
    T = sol.t
    C = sol'
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
            micd[i] = make_MicData(micd[i].MID, micd[i].PID, micd[i].ImT, ti)
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
    if ims == 0
        traj[1] = C[:, :]
    end
    # Delete extinct species
    ms = deleteat!(ms, dls)
    # Now loop over for every immigration attempt
    for i in 1:ims
        if i % 50 == 0
            println("Immigration attempt $i initiated")
            flush(stdout)
        end
        # Find how many immigrants there are
        nI = 1 + rand(sd)
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
        if i != ims
            # Choose a random time for next immigration
            ti = rand(td)
            # Save this as the immigration time
            its[i + 1] = ti + its[i]
            # Then use this and previous immigration time to define the time step
            tspan = (its[i], its[i + 1])
        else
            # At last step just integrate for five times the average time, so that dynamics settle
            tf = 5 * mT + its[i]
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
        C = sol'
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
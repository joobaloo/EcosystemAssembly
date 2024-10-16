# A script to run the dynamics for the full model.
export full_simulate, sing_pop, doub_pop, θ, θ_smooth, qs

# These are temporarily being output to aid with testing
export γs, λs, Eα

# function to find the thermodynamic term θ, for the case of 1 to 1 stoichiometry
function θ(S::Float64, P::Float64, T::Float64, η::Float64, ΔG0::Float64)
    # Catch perverse cases that sometimes arise
    if S <= 0.0
        θs = 1.0
    elseif P <= 0.0
        θs = 0.0
    else
        θs = Q(S, P) / Keq(T, η, ΔG0)
    end
    # θ can be greater than 1, this does not have any impact as q cannot be negative
    return (θs)
end

# version of θ function, that smooths output by setting values > 1 to 1
function θ_smooth(S::Float64, P::Float64, T::Float64, η::Float64, ΔG0::Float64)
    # Catch perverse cases that sometimes arise
    if S <= 0.0
        θs = 1.0
    elseif P <= 0.0
        θs = 0.0
    else
        θs = Q(S, P) / Keq(T, η, ΔG0)
    end
    # In this case don't want to return θ values greater than 1
    return (min(θs, 1.0))
end

# Extra function to find rate from thermodynamic inhibition, just useful for plotting
function qs(ps::Microbe, S::Float64, P::Float64, E::Float64, θs::Float64)
    q = ps.kc[1] * E * S * (1 - θs) / (ps.KS[1] + S * (1 + ps.kr[1] * θs))
    return (max(q, 0.0))
end

# function to find the rate of substrate consumption by a particular reaction
function qs(S::Float64,
        P::Float64,
        E::Float64,
        i::Int64,
        ps::Microbe,
        T::Float64,
        r::Reaction)
    # To speed things I don't have a check here to ensure that r.ID matches ps.Reac[i]
    # This is something to check if I start getting errors
    θs = θ(S, P, T, ps.η[i], r.ΔG0)
    q = ps.kc[i] * E * S * (1 - θs) / (ps.KS[i] + S * (1 + ps.kr[i] * θs))
    # Ensure that negative value cannot be returned
    return (max(q, 0.0))
end

# function to calculate the amount of a particular enzyme a strain has
function Eα(ϕR::Float64, ps::Microbe, i::Int64)
    E = ps.MC * (1 - ϕR - ps.ϕH) * ps.ϕP[i] / (ps.n[2 + i])
    return (E)
end

# function to find (energy use dependent) elongation rate γ
function γs(a::Float64, ps::Microbe)
    γ = ps.γm * a / (a + ps.Kγ)
    return (γ)
end

# function to find the growth rate λ
function λs(a::Float64, ϕR::Float64, ps::Microbe)
    # Find elongation rate
    γ = γs(a, ps)
    λ = (γ * ϕR * ps.Pb) / ps.n[1]
    return (λ)
end

# function to find ϕR based on the energy concentration
function ϕ_R(a::Float64, ps::Microbe)
    ϕ = ps.ω * (1 - ps.ϕH) * a / (ps.KΩ + a)
    return (ϕ)
end

# function to implement the consumer resource dynamics
function full_dynamics!(dx::Array{Float64, 1},
        x::Array{Float64, 1},
        ms::Array{Microbe, 1},
        ps::TOParameters,
        rate::Array{Float64, 2},
        t::Float64)
    # loop over the reactions to find reaction rate for each reaction for each strain
    for j in 1:(ps.O)
        # Find substrate and product for this reaction
        for i in 1:length(ms)
            # Check if microbe i performs reaction j
            if j ∈ ms[i].Reacs
                # Find index of this reaction in microbe
                k = findfirst(x -> x == j, ms[i].Reacs)
                # Find amount of enzyme E
                E = Eα(x[2 * length(ms) + ps.M + i], ms[i], k)
                # Then finally calculate reaction rate
                rate[i, j] = qs(x[length(ms) + ps.reacs[j].Rct],
                    x[length(ms) + ps.reacs[j].Prd],
                    E,
                    k,
                    ms[i],
                    ps.T,
                    ps.reacs[ms[i].Reacs[k]])
            else
                rate[i, j] = 0.0
            end
        end
    end
    # Now want to use the rate matrix in the consumer dynamics
    for i in 1:length(ms)
        # Check if strain is effectively extinct
        if x[i] <= 1e-5
            # If so x should be set to zero and should not change from that
            dx[i] = 0.0
            x[i] = 0.0
            # In this case the energy concentration should also be fixed to zero
            dx[length(ms) + ps.M + i] = 0.0
            x[length(ms) + ps.M + i] = 0.0
            # Corresponding proteome fraction also shouldn't shift
            dx[2 * length(ms) + ps.M + i] = 0.0
        else
            # find growth rate for strains that aren't extinct
            λ = λs(x[length(ms) + ps.M + i], x[2 * length(ms) + ps.M + i], ms[i])
            # (growth rate - death rate)*population
            dx[i] = (λ - ms[i].d) * x[i]
            # Now find optimal ribosome fraction
            ϕR = ϕ_R(x[length(ms) + ps.M + i], ms[i])
            # This introduces a time delay
            τ = ms[i].fd / λ
            # Then update actual ribosome fraction
            dx[2 * length(ms) + ps.M + i] = (ϕR - x[2 * length(ms) + ps.M + i]) / τ
            # Energy intake is zero
            J = 0
            # Loop over all reactions to find energy gained by them
            for j in 1:(ms[i].R)
                J += ms[i].η[j] * rate[i, ms[i].Reacs[j]]
            end
            # Add energy intake and subtract translation and dilution from the energy concentration
            dx[length(ms) + ps.M + i] = J -
                                        (ms[i].MC * ms[i].χl + x[length(ms) + ps.M + i]) * λ
        end
    end
    # Do basic resource dynamics
    for i in (length(ms) + 1):(length(ms) + ps.M)
        # fist add external supply of resource and decay
        dx[i] = ps.κ[i - length(ms)] - ps.δ[i - length(ms)] * x[i]
    end
    # Then loop over microbes
    for i in 1:length(ms)
        # Loop over reactions for specific microbe
        for j in 1:(ms[i].R)
            # Increase the product
            dx[length(ms) + ps.reacs[ms[i].Reacs[j]].Prd] += rate[i, ms[i].Reacs[j]] *
                                                             x[i] / NA
            # and decrease the reactant
            dx[length(ms) + ps.reacs[ms[i].Reacs[j]].Rct] -= rate[i, ms[i].Reacs[j]] *
                                                             x[i] / NA
        end
    end
    # Final step to correct for any concentrations that have dropped below threshold (1e-15)
    for i in (length(ms) + 1):(length(ms) + ps.M)
        if x[i] < 1e-15
            x[i] = 1e-15
            dx[i] = 0.0
        end
    end
    # Any ATP numbers that have gone below 0.33 should be removed
    for i in (length(ms) + ps.M + 1):(2 * length(ms) + ps.M)
        if x[i] < 0.33
            x[i] = 0.0
            dx[i] = 0.0
        end
    end
    return (dx)
end

# Simulation code to run one instance of the simulation with a user defined starting
# condition. ps is parameter set, Tmax is the time to integrate to, pop, conc, as and ϕs
# are the initial conditions, mpl is a pool of microbes, mT is mean immigration time,
# ims is the number of immigrations, λIm controls rate of additional immigrants
function full_simulate(ps::TOParameters,
        pop::Float64,
        conc::Float64,
        as::Float64,
        ϕs::Float64,
        mpl::Array{Microbe, 1},
        Ni::Int64,
        mT::Float64,
        ims::Int64,
        λIm::Float64)
    # Preallocate immigration times
    its = zeros(ims)
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

# function to test for single population growth
function sing_pop(ps::TOParameters,
        pop::Float64,
        conc::Float64,
        as::Float64,
        ϕs::Float64,
        mic::Microbe,
        Tmax::Float64)
    # Preallocate memory
    rate = zeros(1, ps.O)
    # Now substitute preallocated memory in
    dyns!(dx, x, ms, t) = full_dynamics!(dx, x, ms, ps, rate, t)
    # Find time span for this step
    tspan = (0, Tmax)
    # Make appropriate initial condition
    concs = [conc, 0.0]
    x0 = [pop; concs; as; ϕs]
    # Then setup and solve the problem
    prob = ODEProblem(dyns!, x0, tspan, [mic])
    sol = DifferentialEquations.solve(prob)
    return (sol', sol.t)
end

# function to test for competition between two populations
function doub_pop(ps::TOParameters,
        pop::Float64,
        conc::Float64,
        as::Float64,
        ϕs::Float64,
        mics::Array{Microbe, 1},
        Tmax::Float64)
    # Check correct number of microbes has been provided
    if length(mics) != 2
        error("must provide two microbes")
    end
    # Preallocate memory
    rate = zeros(2, ps.O)
    # Now substitute preallocated memory in
    dyns!(dx, x, ms, t) = full_dynamics!(dx, x, ms, ps, rate, t)
    # Find time span for this step
    tspan = (0, Tmax)
    # Make appropriate initial condition
    pops = pop * ones(2)
    concs = conc * ones(ps.M)
    ass = as * ones(2)
    ϕss = ϕs * ones(2)
    x0 = [pops; concs; ass; ϕss]
    # Then setup and solve the problem
    prob = ODEProblem(dyns!, x0, tspan, mics)
    sol = DifferentialEquations.solve(prob)
    return (sol', sol.t)
end

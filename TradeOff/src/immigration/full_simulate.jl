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
                       Ni::Int64, # initial number of microbe species
                       mT::Float64, # mean time between immigration events
                       ims::Int64, # number of immigration events
                       λIm::Float64, #number of invading species same as nI kinda
                       sim_length::Float64)  # duration of simulation
                       # Preallocate immigration times
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

    # Time between each immigration event
    imm_step = sim_length/ims
    # Number of microbes in each immigration event 
    imm_mag = λIm

    # If sim_length is about 2 years step each month?
    step = sim_length/24
    steps_vector = 0:step:sim_length

    # initial step
    tspan = (0, step)
    # solve initial problem
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
    # Delete extinct species
    ms = deleteat!(ms, dls)

    for i in 1:steps_vector
        if i % 24 == 0
            println("Month $i initiated")
            flush(stdout)
        end
        nI = imm_mag
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
        
        

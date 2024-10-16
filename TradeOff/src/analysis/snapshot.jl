# Script to calculate features of the trajectories at specific snapshot times
using TradeOff
using JLD

# Function to calculate relevant values at particular time snap shots
function snp_shot()
    # Check that sufficient arguments have been provided
    if length(ARGS) < 3
        error("insufficient inputs provided")
    end
    # Preallocate the variables I want to extract from the input
    rps = 0
    ims = 0
    sim_type = 0
    # Check that all arguments can be converted to integers
    try
        rps = parse(Int64, ARGS[1])
        ims = parse(Int64, ARGS[2])
        sim_type = parse(Int64, ARGS[3])
    catch e
        error("need to provide 3 integers")
    end
    println("Compiled")
    flush(stdout)
    # Load in hardcoded simulation parameters
    Np, Nt, M, d, μrange = sim_paras(sim_type)
    # Define data directory
    data_dir = joinpath(
        pwd(), "Output", "$(Np)Pools$(M)Metabolites$(Nt)Speciesd=$(d)u=$(μrange)")
    # Read in parameter file
    parameter_file = joinpath(data_dir, "Paras$(ims)Ims.jld")
    if ~isfile(parameter_file)
        error("$(ims) immigrations is missing a parameter file")
    end
    # Load parameters
    ps = load(parameter_file, "ps")
    # Preallocate vector to store final times
    Tfs = zeros(rps)
    # Loop over number of repeats to find maximum times
    for i in 1:rps
        # Load in relevant output file
        output_file = joinpath(data_dir, "Run$(i)Data$(ims)Ims.jld")
        if ~isfile(output_file)
            error("$(ims) immigrations run $(i) is missing an output file")
        end
        # Load in microbe data, and immigration times
        T = load(output_file, "T")
        # Store final T value
        Tfs[i] = T[end]
    end
    # Then find maximum final time
    Tmax = maximum(Tfs)
    # Number of steps to calculate stats for
    NumS = 1000
    # Define snap shot times based on this maximum time
    snps = collect(range(0.0, Tmax, length = NumS))
    # Also save the time step for later use
    t_step = snps[2] - snps[1]
    # Preallocate data to save
    ns = zeros(length(snps) - 1, rps)
    gs = zeros(length(snps) - 1, rps)
    stb = zeros(length(snps) - 1, rps)
    inc = zeros(length(snps) - 1, rps)
    dec = zeros(length(snps) - 1, rps)
    st_r = zeros(length(snps) - 1)
    # List of pools already loaded in
    pls = []
    # Array of array to store pools
    pools = Array{Array{Microbe, 1}, 1}(undef, 1)
    # Counter for number of reactions
    NoR = 0
    # Loop over number of repeats
    for i in 1:rps
        # Load in relevant output file
        output_file = joinpath(data_dir, "Run$(i)Data$(ims)Ims.jld")
        if ~isfile(output_file)
            error("$(ims) immigrations run $(i) is missing an output file")
        end
        # Load in microbe data, and immigration times
        T = load(output_file, "T")
        traj = load(output_file, "traj")
        micd = load(output_file, "micd")
        its = load(output_file, "its")
        # Use to construct full trajectory C
        C = merge_data(ps, traj, T, micd, its)
        # Find and save initial population value for this run
        Ni = C[1, 1]
        # Preallocate vector of microbes
        ms = Array{Microbe, 1}(undef, length(micd))
        # Loop over and find each one
        for j in eachindex(micd)
            # check for case where pool hasn't already been loaded in
            if micd[j].PID ∉ pls
                # Add new pool ID in
                pls = cat(pls, micd[j].PID, dims = 1)
                # Find name of pool
                file = joinpath(pwd(), "Pools",
                    "ID=$(micd[j].PID)N=$(Nt)M=$(ps.M)d=$(d)u=$(μrange).jld")
                # Check if this is the first pool
                if length(pls) == 1
                    # If so save the pool
                    pools[1] = load(file, "mics")
                    # Find number of reactions based on this
                    NoR = maximum(pools[1] .↦ :R)
                else
                    # Otherwise just cat it on existing vector
                    pool = load(file, "mics")
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
        # loop over time snapshots
        for j in 1:(length(snps) - 1)
            # Find indices of the time point before and after the snapshot point
            ind1 = findfirst(x -> x >= snps[j], T)
            ind2 = findfirst(x -> x >= snps[j + 1], T)
            # Check that the simulation is still running at this point
            if ~isnothing(ind1)
                st_r[j] += 1
                # Check for new species that entered the system in this time window
                migs = findall(x -> snps[j] <= x < snps[j + 1], micd .↦ :ImT)
                # Count number of new immigrants
                ns[j, i] = length(migs)
                # Also need to calculate the number that grow (over the snapshot period)
                for k in eachindex(migs)
                    # Check that population has increased from the initial value
                    if C[ind2, migs[k]] > Ni
                        gs[j, i] += 1
                    end
                end
                # Find strains that either exist at the start of the snapshot, or immigrate within it
                survs = findall(((C[ind1, 1:length(micd)] .!== 0.0) .&
                                 (.~isnan.(C[ind1, 1:length(micd)]))) .|
                                (snps[j] .<= (micd .↦ :ImT) .< snps[j + 1]))
                # Calculate change necessary to be considered to be growing
                pchng = exp(0.01 * d * t_step) - 1.00
                # Adjust for the fact that ind2 can be nothing
                if ~isnothing(ind2)
                    indf = ind2
                else
                    indf = lastindex(T)
                end
                # Loop over survivors
                for k in eachindex(survs)
                    # catch species that are extinct by the end of the time window
                    if isnan(C[indf, survs[k]]) || C[indf, survs[k]] == 0.0
                        dec[j, i] += 1
                        # Are they a new immigrant
                    elseif C[ind1, survs[k]] == 0.0 || isnan(C[ind1, survs[k]])
                        # If so do they increase in value
                        if C[ind1, survs[k]] >= Ni
                            inc[j, i] += 1
                            # Other wise they can be treated as decreasing
                        else
                            dec[j, i] += 1
                        end
                        # then ones that have declined sufficiently
                    elseif (C[indf, survs[k]] - C[ind1, survs[k]]) /
                           abs(C[ind1, survs[k]]) < -pchng
                        dec[j, i] += 1
                        # increased sufficiently
                    elseif (C[indf, survs[k]] - C[ind1, survs[k]]) /
                           abs(C[ind1, survs[k]]) > pchng
                        inc[j, i] += 1
                        # If none of the above true assign to stable
                    else
                        stb[j, i] += 1
                    end
                end
            end
        end
        println("Run $i analysed")
        flush(stdout)
    end
    # Now just save the relevant data
    jldopen(joinpath(data_dir, "SnapData$(ims)Ims.jld"), "w") do file
        # Save times of snapshots
        write(file, "times", snps)
        # Save whatever I generate here
        write(file, "ns", ns)
        write(file, "gs", gs)
        write(file, "stb", stb)
        write(file, "inc", inc)
        write(file, "dec", dec)
        write(file, "st_r", st_r)
    end
    return (nothing)
end

@time snp_shot()

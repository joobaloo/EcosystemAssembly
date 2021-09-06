# Script to find how variables change over time, which then saves them
using TradeOff
using JLD

function v_over_t()
    # Check that sufficent arguments have been provided
    if length(ARGS) < 3
        error("insufficent inputs provided")
    end
    # Preallocate the variables I want to extract from the input
    rps = 0
    ims = 0
    sim_type = 0
    # Check that all arguments can be converted to integers
    try
        rps = parse(Int64,ARGS[1])
        ims = parse(Int64,ARGS[2])
        sim_type = parse(Int64,ARGS[3])
    catch e
            error("need to provide 3 integers")
    end
    println("Compiled")
    flush(stdout)
    # Load in hardcoded simulation parameters
    Np, Nt, M, d, μrange = sim_paras(sim_type)
    # Read in parameter file
    pfile = "Output/$(Np)Pools$(M)Metabolites$(Nt)Speciesd=$(d)u=$(μrange)/Paras$(ims)Ims.jld"
    if ~isfile(pfile)
        error("$(ims) immigrations is missing a parameter file")
    end
    # Load parameters
    ps = load(pfile,"ps")
    # List of pools already loaded in
    pls = []
    # Array of array to store pools
    pools = Array{Array{Microbe,1},1}(undef,1)
    # Counter for number of reactions
    NoR = 0
    # Loop over number of repeats
    for i = 1:rps
        # Load in relevant output file
        ofile = "Output/$(Np)Pools$(M)Metabolites$(Nt)Speciesd=$(d)u=$(μrange)/Run$(i)Data$(ims)Ims.jld"
        if ~isfile(ofile)
            error("$(ims) immigrations run $(rN) is missing an output file")
        end
        # Load in microbe data, and immigration times
        T = load(ofile,"T")
        traj = load(ofile,"traj")
        micd = load(ofile,"micd")
        its = load(ofile,"its")
        # Use to construct full trajectory C
        C = merge_data(ps,traj,T,micd,its)
        # Preallocate vector of microbes
        ms = Array{Microbe,1}(undef,length(micd))
        # Loop over and find each one
        for j = 1:length(micd)
            # check for case where pool hasn't already been loaded in
            if micd[j].PID ∉ pls
                # Add new pool ID in
                pls = cat(pls,micd[j].PID,dims=1)
                # Find name of pool
                file = "Pools/ID=$(micd[j].PID)N=$(Nt)M=$(ps.M)d=$(d)u=$(μrange).jld"
                # Check if this is the first pool
                if length(pls) == 1
                    # If so save the pool
                    pools[1] = load(file,"mics")
                    # Find number of reactions based on this
                    NoR = maximum(pools[1].↦:R)
                else
                    # Otherwise just cat it on existing vector
                    pool = load(file,"mics")
                    pools = cat(pools,pool,dims=1)
                    # Find maximum number of reactions for this pool
                    NoRt = maximum(pools[1].↦:R)
                    # Save if higher than old number of reactions
                    NoR = max(NoR,NoRt)
                end
            end
            # Find correct pool to read from
            ind = findfirst(x->x==micd[j].PID,pls)
            # Use this index to find and save the correct microbe
            ms[j] = (pools[ind])[micd[j].MID]
        end
        # Preallocate containers to store number of survivors with time
        svt = Array{Int64,1}(undef,length(T))
        tsvt = Array{Int64,1}(undef,length(T))
        sbs = Array{Int64,1}(undef,length(T))
        Rs = Array{Int64,2}(undef,NoR,length(T))
        via_R = Array{Int64,2}(undef,NoR,length(T))
        ηs = zeros(length(T))
        kcs = zeros(length(T))
        KSs = zeros(length(T))
        krs = zeros(length(T))
        av_steps = zeros(length(T))
        via_η = zeros(length(T))
        ωs = zeros(length(T))
        via_ω = zeros(length(T))
        fr_ΔG = zeros(length(T))
        via_a = zeros(length(T))
        via_ϕR = zeros(length(T))
        η1 = zeros(length(T))
        η2 = zeros(length(T))
        fr_ΔG1 = zeros(length(T))
        fr_ΔG2 = zeros(length(T))
        ηs_R = zeros(NoR,length(T))
        ωs_R = zeros(NoR,length(T))
        kc_R = zeros(NoR,length(T))
        KS_R = zeros(NoR,length(T))
        kr_R = zeros(NoR,length(T))
        # Save total number of strains
        numS = length(micd)
        # Make vector of indices
        a_i = collect((numS+ps.M+1):(2*numS+ps.M))
        ϕ_i = collect((2*numS+ps.M+1):(3*numS+ps.M))
        # Loop over all time points
        for j = 1:length(T)
            # Find indices of surviving strains
            inds = findall(x->x>1e-5,C[j,1:numS])
            # Save number of surviving strains at each time point
            svt[j] = length(inds)
            # Find indices of "viable" strains
            vinds = findall(x->x>1e5,C[j,1:numS])
            # Save number of "viable" strains
            tsvt[j] = length(vinds)
            # Then also number of substrates (skipping final waste product)
            sbs[j] = count(x->x>1e-12,C[j,(numS+1):(numS+ps.M-1)])
            # Calculate average energy concentration (a)
            if length(a_i[vinds]) > 0
                via_a[j] = sum(C[j,a_i[vinds]])/length(a_i[vinds])
            end
            # Do the same for the ribosome fraction
            if length(ϕ_i[vinds]) > 0
                via_ϕR[j] = sum(C[j,ϕ_i[vinds]])/length(ϕ_i[vinds])
            end
            # Loop over number of reactions
            for k = 1:NoR
                # Count number of strains with reaction for each case
                Rs[k,j] = count(x->x==k,ms[inds].↦:R)
                via_R[k,j] = count(x->x==k,ms[vinds].↦:R)
            end
            # Find (weighted) total eta value, and ω value, and kinetic parameters
            for k = 1:length(inds)
                ηs[j] += sum(ms[inds[k]].η.*ms[inds[k]].ϕP)
                ωs[j] += ms[inds[k]].ω
            end
            # Average over number of strains
            if svt[j] > 0
                ηs[j] /= svt[j]
                ωs[j] /= svt[j]
            end
            # Set up counters for the number of strains with 1 and 2 gap reactions, respectively
            c1 = 0
            c2 = 0
            # Find (weighted) total eta value for viable strains
            for k = 1:length(vinds)
                via_η[j] += sum(ms[vinds[k]].η.*ms[vinds[k]].ϕP)
                via_ω[j] += ms[vinds[k]].ω
                kcs[j] += sum(ms[vinds[k]].kc.*ms[vinds[k]].ϕP)
                KSs[j] += sum(ms[vinds[k]].KS.*ms[vinds[k]].ϕP)
                krs[j] += sum(ms[vinds[k]].kr.*ms[vinds[k]].ϕP)
                # Bools to store whether strain has 1 gap and 2 gap reaction, respectively
                pres1 = false
                pres2 = false
                # Loop over all reactions this strain has
                for l = 1:ms[vinds[k]].R
                    # Find reaction number
                    Rn = ms[vinds[k]].Reacs[l]
                    # Find relevant reaction
                    r = ps.reacs[Rn]
                    # Then calculate frac transduced
                    fr_ΔG[j] += ms[vinds[k]].η[l].*ms[vinds[k]].ϕP[l]*ΔGATP/(-r.ΔG0)
                    # Find step size
                    s_size = (r.Prd - r.Rct)
                    # weight this step size to 1 and add to total
                    av_steps[j] += (s_size)*ms[vinds[k]].ϕP[l]
                    # use step size to choose which eta value to add to
                    if s_size == 1
                        pres1 = true
                        η1[j] += ms[vinds[k]].η[l].*ms[vinds[k]].ϕP[l]
                        fr_ΔG1[j] += ms[vinds[k]].η[l].*ms[vinds[k]].ϕP[l]*ΔGATP/(-r.ΔG0)
                    elseif s_size == 2
                        pres2 = true
                        η2[j] += ms[vinds[k]].η[l].*ms[vinds[k]].ϕP[l]
                        fr_ΔG2[j] += ms[vinds[k]].η[l].*ms[vinds[k]].ϕP[l]*ΔGATP/(-r.ΔG0)
                    end
                end
                # If they are present them increment the counters
                if pres1 == true
                    c1 += 1
                end
                if pres2 == true
                    c2 += 1
                end
            end
            # Average over number of strains
            if tsvt[j] > 0
                via_η[j] /= tsvt[j]
                via_ω[j] /= tsvt[j]
                kcs[j] /= tsvt[j]
                KSs[j] /= tsvt[j]
                krs[j] /= tsvt[j]
                av_steps[j] /= tsvt[j]
                fr_ΔG[j] /= tsvt[j]
                # Divide by number of strains that possess reactions
                if c1 > 0
                    η1[j] /= c1
                    fr_ΔG1[j] /= c1
                end
                if c2 > 0
                    η2[j] /= c2
                    fr_ΔG2[j] /= c2
                end
            end
            # Break down eta and omega value by R
            for k = 1:length(vinds)
                # Find relevant reaction number
                l = ms[vinds[k]].R
                # Add contribution to relevant total
                ηs_R[l,j] += sum(ms[vinds[k]].η.*ms[vinds[k]].ϕP)
                ωs_R[l,j] += ms[vinds[k]].ω
                kc_R[l,j] += sum(ms[vinds[k]].kc.*ms[vinds[k]].ϕP)
                KS_R[l,j] += sum(ms[vinds[k]].KS.*ms[vinds[k]].ϕP)
                kr_R[l,j] += sum(ms[vinds[k]].kr.*ms[vinds[k]].ϕP)
            end
            # Now weight by number of strains with each type of reaction
            for k = 1:NoR
                if via_R[k,j] > 0
                    ηs_R[k,j] /= via_R[k,j]
                    ωs_R[k,j] /= via_R[k,j]
                    kc_R[k,j] /= via_R[k,j]
                    KS_R[k,j] /= via_R[k,j]
                    kr_R[k,j] /= via_R[k,j]
                end
            end
        end
        # Now just save the relevant data
        jldopen("Output/$(Np)Pools$(M)Metabolites$(Nt)Speciesd=$(d)u=$(μrange)/AvRun$(i)Data$(ims)Ims.jld","w") do file
            # Save full timecourse
            write(file,"T",T)
            # Save reaction data
            write(file,"Rs",Rs)
            write(file,"via_R",via_R)
            write(file,"ηs_R",ηs_R)
            write(file,"ωs_R",ωs_R)
            write(file,"kc_R",kc_R)
            write(file,"KS_R",KS_R)
            write(file,"kr_R",kr_R)
            # Save the other quantities
            write(file,"svt",svt)
            write(file,"tsvt",tsvt)
            write(file,"sbs",sbs)
            write(file,"ηs",ηs)
            write(file,"via_η",via_η)
            write(file,"kcs",kcs)
            write(file,"KSs",KSs)
            write(file,"krs",krs)
            write(file,"av_steps",av_steps)
            write(file,"ωs",ωs)
            write(file,"via_ω",via_ω)
            write(file,"fr_ΔG",fr_ΔG)
            write(file,"via_a",via_a)
            write(file,"via_ϕR",via_ϕR)
            write(file,"η1",η1)
            write(file,"η2",η2)
            write(file,"fr_ΔG1",fr_ΔG1)
            write(file,"fr_ΔG2",fr_ΔG2)
            # Finally save final time to help with benchmarking
            write(file,"Tf",T[end])
        end
        println("Run $i analysed")
        flush(stdout)
    end
    return(nothing)
end


@time v_over_t()

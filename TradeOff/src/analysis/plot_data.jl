# A script to read in and analyse model output data
using TradeOff
using JLD
using Plots
import PyPlot

# function to plot the trajectories
function plot_traj()
    # Check that sufficient arguments have been provided
    if length(ARGS) < 3
        error("insufficient inputs provided")
    end
    # Preallocate the variables I want to extract from the input
    rN = 0
    ims = 0
    sim_type = 0
    # Check that all arguments can be converted to integers
    try
        rN = parse(Int64, ARGS[1])
        ims = parse(Int64, ARGS[2])
        sim_type = parse(Int64, ARGS[3])
    catch e
        error("need to provide 3 integers")
    end
    println("Compiled")
    # Extract other simulation parameters from the function
    Np, Nt, M, d, μrange = sim_paras(sim_type)
    # Token to insert into filenames
    tk = ""
    # Overwritten for no immigration case
    if sim_type == 5
        tk = "NoImm"
    end
    # Define data directory to read from
    data_dir = joinpath(
        pwd(), "Output", "$(tk)$(Np)Pools$(M)Metabolites$(Nt)Speciesd=$(d)u=$(μrange)")
    # Read in appropriate files
    parameter_file = joinpath(data_dir, "Paras$(ims)Ims.jld")
    if ~isfile(parameter_file)
        error("$(ims) immigrations run $(rN) is missing a parameter file")
    end
    output_file = joinpath(data_dir, "Run$(rN)Data$(ims)Ims.jld")
    if ~isfile(output_file)
        error("$(ims) immigrations run $(rN) is missing an output file")
    end
    # Read in relevant data
    ps = load(parameter_file, "ps")
    traj = load(output_file, "traj")
    T = load(output_file, "T")
    micd = load(output_file, "micd")
    its = load(output_file, "its")
    println("Data read in")
    # Find C from a function
    C = merge_data(ps, traj, T, micd, its)
    println("Data merged")
    # Define and (if necessary create) the directory for the plots
    plot_dir = joinpath(pwd(), "Output", "$(tk)Plotsd=$(d)u=$(μrange)")
    mkpath(plot_dir)
    # Find total number of strains
    totN = length(micd)
    # Find indices of extinct strains
    ext = isnan.(C[end, 1:totN])
    # Invert to find survivors
    svr = .!ext
    pyplot(dpi = 200)
    # Plot all the populations
    p1 = plot(yaxis = :log10, ylabel = "Population (# cells)", ylims = (1e-5, Inf))
    for i in 1:totN
        # Find and eliminate zeros so that they can be plotted on a log plot
        inds = (C[:, i] .> 0)
        plot!(p1, T[inds], C[inds, i], label = "")
    end
    savefig(p1, joinpath(plot_dir, "all_pops.png"))
    # Plot all the concentrations
    p2 = plot(yaxis = :log10, ylabel = "Concentration")#,ylims=(1e-15,Inf))
    for i in 1:(ps.M)
        # Find and eliminate zeros so that they can be plotted on a log plot
        inds = (C[:, totN + i] .> 0)
        plot!(p2, T[inds], C[inds, totN + i], label = "")
    end
    savefig(p2, joinpath(plot_dir, "all_concs.png"))
    # Plot all the energy concentrations
    p3 = plot(ylabel = "Energy Concentration")
    for i in 1:totN
        plot!(p3, T, C[:, totN + ps.M + i], label = "")
    end
    savefig(p3, joinpath(plot_dir, "all_as.png"))
    # Plot all the ribosome fractions
    p4 = plot(ylabel = "Ribosome fraction")
    for i in 1:totN
        plot!(p4, T, C[:, 2 * totN + ps.M + i], label = "")
    end
    savefig(p4, joinpath(plot_dir, "all_fracs.png"))
    # Plot populations that survive to the end
    p1 = plot(yaxis = :log10, ylabel = "Population (# cells)", ylims = (1e-5, Inf))
    for i in 1:totN
        if svr[i] == true
            # Find and eliminate zeros so that they can be plotted on a log plot
            inds = (C[:, i] .> 0)
            plot!(p1, T[inds], C[inds, i], label = "")
        end
    end
    savefig(p1, joinpath(plot_dir, "surv_pops.png"))
    # Plot energy concentrations of populations that survive to the end
    p3 = plot(ylabel = "Energy Concentration")
    for i in 1:totN
        if svr[i] == true
            plot!(p3, T, C[:, totN + ps.M + i], label = "")
        end
    end
    savefig(p3, joinpath(plot_dir, "surv_as.png"))
    # Plot ribosome fractions of populations that survive to the end
    p4 = plot(ylabel = "Ribosome fraction")
    for i in 1:totN
        if svr[i] == true
            plot!(p4, T, C[:, 2 * totN + ps.M + i], label = "")
        end
    end
    savefig(p4, joinpath(plot_dir, "surv_fracs.png"))
    return (nothing)
end

# function to plot the averages for just one run
function plot_run_averages()
    # Check that sufficient arguments have been provided
    if length(ARGS) < 3
        error("insufficient inputs provided")
    end
    # Preallocate the variables I want to extract from the input
    ims = 0
    rN = 0
    sim_type = 0
    # Check that all arguments can be converted to integers
    try
        rN = parse(Int64, ARGS[1])
        ims = parse(Int64, ARGS[2])
        sim_type = parse(Int64, ARGS[3])
    catch e
        error("need to provide three integers")
    end
    println("Compiled")
    # Load in hardcoded simulation parameters
    Np, Nt, M, d, μrange = sim_paras(sim_type)
    output_dir = joinpath(pwd(), "Output")
    averages_file = joinpath(
        output_dir, "$(Np)Pools$(M)Metabolites$(Nt)Speciesd=$(d)u=$(μrange)",
        "AvRun$(rN)Data$(ims)Ims.jld")
    if ~isfile(averages_file)
        error("$(ims) immigrations run $(rN) is missing an output file")
    end
    # Read in relevant data
    T = load(averages_file, "T")
    svt = load(averages_file, "surviving_species")
    tsvt = load(averages_file, "viable_species")
    ηs = load(averages_file, "average_η")
    # Setup plotting
    pyplot(dpi = 200)
    # Plot this data
    plot(T, svt, label = "", xlabel = "Time (s)", ylabel = "Number of survivors")
    savefig(joinpath(output_dir, "SvTime.png"))
    plot(T, ηs, label = "", xlabel = "Time (s)", ylabel = "eta")
    savefig(joinpath(output_dir, "EtaTime.png"))
    plot(T, tsvt, label = "", xlabel = "Time (s)", ylabel = "Number of viable strains")
    savefig(joinpath(output_dir, "ViaTime.png"))
    return (nothing)
end

@time plot_run_averages()
@time plot_traj()

# Script to plot figure 4
using TradeOff
using Plots
using JLD
using ColorSchemes
using LaTeXStrings
using Plots.PlotMeasures

function sub_prob(M::Int64, R::Int64, subs::Float64)
    # Check that final waste product hasn't been generated
    if subs >= M - 1
        no_sub = M - 1

    elseif subs <= 0.0 # Also check that it hasn't somehow gone negative
        no_sub = 0.0
    else
        no_sub = subs
    end
    # Calculate probability
    P = 1 - (1 - (no_sub) / (M - 1))^R
    return (P)
end

function figure4(ims::Int64, sim_type::Int64, sim_type2::Int64)
    println("Compiled")
    # Extract other simulation parameters from the function
    Np, Nt, M, d, μrange = sim_paras(sim_type)
    # Find file name to load in
    data_dir = joinpath(
        pwd(), "Output", "$(Np)Pools$(M)Metabolites$(Nt)Speciesd=$(d)u=$(μrange)")
    stats_file = joinpath(data_dir, "RunStats$(ims)Ims.jld")
    # Check it actually exists
    if ~isfile(stats_file)
        error("missing stats file for $(ims) immigrations simulations")
    end
    # Now load out the times, and number of trajectories
    times = load(stats_file, "times")
    no_sims = load(stats_file, "no_sims")
    no_via = load(stats_file, "no_via")
    no_rs = load(stats_file, "no_rs")
    # Load in averages
    mn_sbs = load(stats_file, "mn_sbs")
    mn_via_R = load(stats_file, "mn_via_R")
    mn_ηs_R = load(stats_file, "mn_ηs_R")
    mn_KS_R = load(stats_file, "mn_KS_R")
    # Load in standard deviations
    sd_sbs = load(stats_file, "sd_sbs")
    sd_via_R = load(stats_file, "sd_via_R")
    sd_ηs_R = load(stats_file, "sd_ηs_R")
    sd_KS_R = load(stats_file, "sd_KS_R")
    # Preallocate standard errors
    se_via_R = zeros(size(sd_via_R))
    se_ηs_R = zeros(size(sd_ηs_R))
    se_KS_R = zeros(size(sd_KS_R))
    # Calculate standard errors from this
    se_sbs = sd_sbs ./ sqrt.(no_sims)
    # Calculation (slightly) different in the viable case
    for i in axes(sd_via_R, 1)
        se_via_R[i, :] = sd_via_R[i, :] ./ sqrt.(no_via)
        se_ηs_R[i, :] = sd_ηs_R[i, :] ./ sqrt.(no_rs[i, :])
        se_KS_R[i, :] = sd_KS_R[i, :] ./ sqrt.(no_rs[i, :])
    end
    # Preallocate probabilities
    mn_Ps = zeros(size(mn_via_R))
    up_Ps = zeros(size(mn_via_R))
    dw_Ps = zeros(size(mn_via_R))
    # Loop over, calculating the probability at each step
    for i in axes(mn_Ps, 1)
        for j in axes(mn_Ps, 2)
            # Calculate mean probability
            mn_Ps[i, j] = sub_prob(M, i, mn_sbs[j])
            # And also upper and lower bound
            up_Ps[i, j] = sub_prob(M, i, mn_sbs[j] + se_sbs[j]) .- mn_Ps[i, j]
            dw_Ps[i, j] = mn_Ps[i, j] .- sub_prob(M, i, mn_sbs[j] - se_sbs[j])
        end
    end
    println("Data read in")
    # Define output directory and if necessary make it
    outdir = joinpath(pwd(), "Output", "Fig4")
    mkpath(outdir)
    # Set default plotting options
    default(dpi = 200)
    # Load in colour scheme
    a = ColorSchemes.sunset.colors
    # Plot basic trade-off first
    p1 = plot(xlabel = "Time (s)",
        ylabel = "Number of species",
        xlim = (-Inf, 5e7),
        title = "High substrate free-energy case",
        legend = :bottomright,
        ylim = (0.0, 8.0))
    plot!(p1, times, mn_via_R[1, :], ribbon = se_via_R[1, :], label = "R=1", color = a[1])
    plot!(p1, times, mn_via_R[3, :], ribbon = se_via_R[3, :], label = "R=3", color = a[2])
    plot!(p1, times, mn_via_R[5, :], ribbon = se_via_R[5, :], label = "R=5", color = a[3])
    plot!(p1, times, mn_via_R[7, :], ribbon = se_via_R[7, :], label = "R=7", color = a[4])
    # Add annotation
    px, py = annpos([0.0; 5e7], [0.0; 8.0], 0.075, 0.05)
    annotate!(p1, px, py, text("A", 17, :black))
    savefig(p1, joinpath(outdir, "AvViaReacsTime.png"))
    # Now do probability plot
    p2 = plot(xlabel = "Time (s)",
        ylabel = "Probability of no usable substrate",
        xlim = (-Inf, 5e7),
        title = "Chance of species finding no usable substrates",
        legend = false)
    plot!(p2,
        times,
        1 .- mn_Ps[1, :],
        ribbon = (up_Ps[1, :], dw_Ps[1, :]),
        label = "R=1",
        color = a[1])
    plot!(p2,
        times,
        1 .- mn_Ps[3, :],
        ribbon = (up_Ps[3, :], dw_Ps[3, :]),
        label = "R=3",
        color = a[2])
    plot!(p2,
        times,
        1 .- mn_Ps[5, :],
        ribbon = (up_Ps[5, :], dw_Ps[5, :]),
        label = "R=5",
        color = a[3])
    plot!(p2,
        times,
        1 .- mn_Ps[7, :],
        ribbon = (up_Ps[7, :], dw_Ps[7, :]),
        label = "R=7",
        color = a[4])
    # Add annotation
    px, py = annpos([0.0; 5e7], [0.0; 1.0375], 0.075, 0.05)
    annotate!(p2, px, py, text("B", 17, :black))
    savefig(p2, joinpath(data_dir, "ProbSubTime.png"))
    # Plot trade-off for η
    p3 = plot(xlabel = "Time (s)",
        ylabel = "Average eta value",
        xlim = (-Inf, 5e7),
        title = "Variation of key parameters",
        legend = false,
        ylim = (0.0, 7.5))
    plot!(p3, times, mn_ηs_R[1, :], ribbon = se_ηs_R[1, :], label = "R=1", color = a[1])
    plot!(p3, times, mn_ηs_R[3, :], ribbon = se_ηs_R[3, :], label = "R=3", color = a[2])
    plot!(p3, times, mn_ηs_R[5, :], ribbon = se_ηs_R[5, :], label = "R=5", color = a[3])
    plot!(p3, times, mn_ηs_R[7, :], ribbon = se_ηs_R[7, :], label = "R=7", color = a[4])
    # Add annotation
    px, py = annpos([0.0; 5e7], [0.0; 7.5], 0.075, 0.05)
    annotate!(p3, px, py, text("C", 17, :black))
    # Define box for inset here
    box = (1, bbox(0.4, 0.15, 0.4, 0.3, :bottom, :left))
    Ks = L"K_S"
    e7 = L"10^7"
    em3 = L"10^{-3}"
    # Plot other trade-off into the inset
    plot!(p3,
        times / 1e7,
        mn_KS_R[1, :] * 1000.0,
        ribbon = se_KS_R[1, :] * 1000.0,
        color = a[1],
        label = "",
        inset_subplots = box,
        subplot = 2)
    plot!(p3,
        times / 1e7,
        mn_KS_R[3, :] * 1000.0,
        ribbon = se_KS_R[3, :] * 1000.0,
        color = a[2],
        label = "",
        subplot = 2)
    plot!(p3,
        times / 1e7,
        mn_KS_R[5, :] * 1000.0,
        ribbon = se_KS_R[5, :] * 1000.0,
        color = a[3],
        label = "",
        subplot = 2)
    plot!(p3,
        times / 1e7,
        mn_KS_R[7, :] * 1000.0,
        ribbon = se_KS_R[7, :] * 1000.0,
        color = a[4],
        label = "",
        subplot = 2)
    plot!(p3,
        xlabel = "Time ($(e7) s)",
        ylabel = "$(Ks) ($em3)",
        xlim = (-Inf, 5.0),
        grid = false,
        subplot = 2)
    savefig(p3, joinpath(data_dir, "AvEtaperReacTime.png"))
    # Extract simulation parameters for the other case from the function
    Np, Nt, M, d, μrange = sim_paras(sim_type2)
    # Find file name to load in
    stats_file_2 = joinpath(data_dir, "RunStats$(ims)Ims.jld")
    # Check it actually exists
    if ~isfile(stats_file_2)
        error("missing stats file for $(ims) immigrations simulations")
    end
    # Now load out the times, and number of trajectories
    times = load(stats_file_2, "times")
    no_via = load(stats_file_2, "no_via")
    # Load in averages
    mn_via_R = load(stats_file_2, "mn_via_R")
    # Load in standard deviations
    sd_via_R = load(stats_file_2, "sd_via_R")
    # Preallocate standard errors
    se_via_R = zeros(size(sd_via_R))
    # Calculation (slightly) different in the viable case
    for i in axes(sd_via_R, 1)
        se_via_R[i, :] = sd_via_R[i, :] ./ sqrt.(no_via)
    end
    p4 = plot(xlabel = "Time (s)",
        ylabel = "Number of species",
        xlim = (-Inf, 5e7),
        title = "Low substrate free-energy case",
        legend = false,
        ylim = (0.0, 5.0))
    plot!(p4, times, mn_via_R[1, :], ribbon = se_via_R[1, :], label = "R=1", color = a[1])
    plot!(p4, times, mn_via_R[3, :], ribbon = se_via_R[3, :], label = "R=3", color = a[2])
    plot!(p4, times, mn_via_R[5, :], ribbon = se_via_R[5, :], label = "R=5", color = a[3])
    plot!(p4, times, mn_via_R[7, :], ribbon = se_via_R[7, :], label = "R=7", color = a[4])
    # Add annotation
    px, py = annpos([0.0; 5e7], [0.0; 5.0], 0.075, 0.05)
    annotate!(p4, px, py, text("D", 17, :black))
    savefig(p4, joinpath(outdir, "LowFreeEnergy.png"))
    # Plot all graphs as a single figure
    pt = plot(p1, p3, p2, p4, layout = 4, size = (1200, 800), margin = 5.0mm)
    savefig(pt, joinpath(outdir, "figure4.png"))
    return (nothing)
end

@time figure4(500, 1, 2)

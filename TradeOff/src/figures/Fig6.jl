# Script to plot figure 6
using TradeOff
using Plots
using JLD
using ColorSchemes
using Plots.PlotMeasures
using LaTeXStrings

function find_label(sim_type::Int64)
    # Assign label based on simulation type
    if sim_type == 1
        lb = "high free-energy"
    else
        lb = "low free-energy"
    end
    return (lb)
end

function figure6(ims::Int64)
    println("Compiled")
    # Set default plotting options
    default(dpi = 200)
    # Load in colour schemes
    a = ColorSchemes.Dark2_4.colors
    # Extract specific 4 colours from a colour scheme
    bt = ColorSchemes.tab10.colors
    b = [bt[10]; bt[9]; bt[7]; bt[5]]
    # Make plot objects
    p1 = plot(xlabel = "Times (s)",
        xlim = (-Inf, 5e7),
        legend = :right,
        ylabel = L"\eta",
        title = "ATP yield",
        ylim = (1.0, 6.0))
    p2 = plot(xlabel = "Times (s)",
        xlim = (-Inf, 5e7),
        ylim = (0.55, 0.75),
        legend = false,
        ylabel = "Fraction of free-energy transduced",
        title = "Average reaction efficiency")
    p3 = plot(xlabel = "Times (s)",
        xlim = (-Inf, 5e7),
        ylim = (2.25, 3.00),
        legend = false,
        ylabel = "Average number of reaction steps",
        title = "Relative frequency of reaction types")
    # Loop over the 4 conditions
    for i in 1:2
        # Extract other simulation parameters from the function
        Np, Nt, M, d, μrange = sim_paras(i)
        data_dir = joinpath(
            pwd(), "Output", "$(Np)Pools$(M)Metabolites$(Nt)Speciesd=$(d)u=$(μrange)")
        # Find file name to load in
        stats_file = joinpath(data_dir, "RunStats$(ims)Ims.jld")
        # Check it actually exists
        if ~isfile(stats_file)
            error("missing stats file for $(ims) immigrations simulations")
        end
        # Now load out the times, and number of trajectories
        times = load(stats_file, "times")
        no_via = load(stats_file, "no_via")
        # Load in averages
        mn_via_η_bw = load(stats_file, "mn_via_η_bw")
        mn_fr_ΔG_bw = load(stats_file, "mn_fr_ΔG_bw")
        mn_av_steps_bw = load(stats_file, "mn_av_steps_bw")
        # Load in standard deviations
        sd_via_η_bw = load(stats_file, "sd_via_η_bw")
        sd_fr_ΔG_bw = load(stats_file, "sd_fr_ΔG_bw")
        sd_av_steps_bw = load(stats_file, "sd_av_steps_bw")
        # Calculate relevant standard errors
        se_via_η_bw = sd_via_η_bw ./ sqrt.(no_via)
        se_fr_ΔG_bw = sd_fr_ΔG_bw ./ sqrt.(no_via)
        se_av_steps_bw = sd_av_steps_bw ./ sqrt.(no_via)
        # Find appropriate label
        lb = find_label(i)
        # Plot the data to the relevant plot objects
        plot!(p1, times, mn_via_η_bw, ribbon = se_via_η_bw, color = a[i], label = lb)
        plot!(p2, times, mn_fr_ΔG_bw, ribbon = se_fr_ΔG_bw, color = a[i])
        plot!(p3, times, mn_av_steps_bw, ribbon = se_av_steps_bw, color = a[i])
    end
    # Load in no immigration data
    Np, Nt, M, d, μrange = sim_paras(5)
    data_dir = joinpath(
        pwd(), "Output", "NoImm$(Np)Pools$(M)Metabolites$(Nt)Speciesd=$(d)u=$(μrange)")
    noim_file = joinpath(data_dir, "RunStats0Ims.jld")
    # Check it actually exists
    if ~isfile(noim_file)
        error("missing stats file for no immigrations simulations")
    end
    times = load(noim_file, "times")
    no_via = load(noim_file, "no_via")
    # Load in η data
    mn_via_η_bw = load(noim_file, "mn_via_η_bw")
    sd_via_η_bw = load(noim_file, "sd_via_η_bw")
    # Calculate relevant standard errors
    se_via_η_bw = sd_via_η_bw ./ sqrt.(no_via)
    # Then plot
    p4 = plot(times,
        mn_via_η_bw,
        ribbon = se_via_η_bw,
        color = a[3],
        label = "high free-energy",
        xlim = (-Inf, 2.5e6),
        ylim = (1.0, 6.0),
        title = "ATP yield without immigration")
    plot!(p4, xlabel = "Times (s)", ylabel = L"\eta")
    # Define output directory and if necessary make it
    outdir = joinpath(pwd(), "Output", "Fig6")
    mkpath(outdir)
    # Save figures to this directory
    savefig(p1, joinpath(outdir, "Eta.png"))
    savefig(p2, joinpath(outdir, "Efficiency.png"))
    savefig(p3, joinpath(outdir, "Steps.png"))
    savefig(p4, joinpath(outdir, "NoImmEta.png"))
    # Add annotations
    px, py = annpos([0.0; 5e7], [1.0; 6.0], 0.075, 0.05)
    annotate!(p1, px, py, text("A", 17, :black))
    px, py = annpos([0.0; 5e7], [0.55; 0.75], 0.075, 0.05)
    annotate!(p2, px, py, text("B", 17, :black))
    px, py = annpos([0.0; 5e7], [2.25; 3.00], 0.075, 0.05)
    annotate!(p3, px, py, text("C", 17, :black))
    px, py = annpos([0.0; 2.5e6], [1.0; 6.0], 0.075, 0.05)
    annotate!(p4, px, py, text("D", 17, :black))
    # Plot all graphs as a single figure
    pt = plot(p1, p2, p3, p4, layout = (2, 2), size = (1200, 800), margin = 5.0mm)
    savefig(pt, joinpath(outdir, "figure6.png"))
    return (nothing)
end

@time figure6(500)

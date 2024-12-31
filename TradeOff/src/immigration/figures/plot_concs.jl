using TradeOff
using JLD2
using Plots
include("../simulation_functions.jl")

function plot_traj()
    for i in 1:20
        # Preallocate the variables I want to extract from the input
        rN = i
        rate = 320
        sim_type = 1
        rl = 1
        ru = 5
        # Extract other simulation parameters from the function
        Np, Nt, M, d, μrange = sim_paras(sim_type)

        # Define data directory to read from
        data_dir = joinpath(pwd(), "Output", "niche_size$(rl)_$(ru)", "$(rate)_a_year_rate")

        # Read in appropriate files
        parameter_file = joinpath(data_dir, "Parameters.jld")
        if ~isfile(parameter_file)
            error("$(rate) immigrations run $(rN) is missing a parameter file")
        end

        output_file = joinpath(data_dir, "Run$(rN)Data.jld")
        if ~isfile(output_file)
            error("$(rate) immigrations run $(rN) is missing an output file")
        end

        # Read in relevant data
        ps = load(parameter_file, "ps")
        traj = load(output_file, "traj")
        T = load(output_file, "T")
        micd = load(output_file, "micd")
        its = load(output_file, "its")
        println("Data read in")
        
        # Convert `its` to Vector{Float64}
        its_vector = collect(its)

        # Find C from a function
        C = imm_merge_data(ps, traj, T, micd, its_vector)
        println("Data merged")

        # Define and (if necessary create) the directory for the plots
        plot_dir = joinpath(pwd(), "Output", "niche_size$(rl)_$(ru)", "1immigrants")
        mkpath(plot_dir)

        # Find total number of strains
        totN = length(micd)

        # Plot all the concentrations
        p2 = plot(yaxis = :log10, ylabel = "Concentration", plot_title = "niche_size$(rl)_$(ru)")
        for i in 1:(ps.M)
            # Find and eliminate zeros so that they can be plotted on a log plot
            inds = (C[:, totN + i] .> 0)
            plot!(p2, T[inds], C[inds, totN + i], label = "")
        end
        savefig(p2, joinpath(plot_dir, "run$(i)niche_size$(rl)_$(ru)all_concs.png"))
end
    return (nothing)
end

@time plot_traj()
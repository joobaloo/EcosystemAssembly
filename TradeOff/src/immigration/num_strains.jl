#write a function that saves the number of strains at each timepoint

using TradeOff
using Plots
using JLD
include("../immigration/simulation_functions.jl")


function calc_num_strains()
    
    #num_events_vector = (10, 100, 200, 300, 400, 500)
    num_events_vector = (10, 100, 200)
    # Define an empty array to store vectors
    num_strains_array = Vector{Vector{Int}}()

    for i in 1:length(num_events_vector)
        rps = 1
        num_immigrations = i
        num_immigrants = 10

        # Read in appropriate files
        pfile = "Output/$(num_events_vector[i])events_$(num_immigrants)immigrants/Parameters.jld"
        println("Output/$(num_events_vector[i])events_$(num_immigrants)immigrants/Parameters.jld")
        if ~isfile(pfile)
            error("$(num_immigrants) immigrations run $(rps) is missing a parameter file")
        end
        
        ofile = "Output/$(num_events_vector[i])events_$(num_immigrants)immigrants/Run$(rps)Data.jld"
        if ~isfile(ofile)
            error("$(num_events_vector[i]) immigration events with (num_immigrants) immigrants run $(rps) is missing an output file")
        end

        # Read in relevant data
        ps = load(pfile, "ps")
        traj = load(ofile, "traj")
        T = load(ofile, "T")
        micd = load(ofile, "micd")
        its = load(ofile, "its")
        println("Data read in")

        # Convert `its` to Vector{Float64}
        its_vector = collect(its)

        # Find C from a function
        C = imm_merge_data(ps, traj, T, micd, its_vector)
        println("Data merged")
        println("calculating species richness for $(num_events_vector[i]) events")

        current_num_strains = sum(C .> 0, dims=2)[:, 1]
        
        # Add the generated vector to the array
        push!(num_strains_array, current_num_strains)
    end


    # Specify the file path
    jld_file_path = "Output/Imm_plots/number_of_strains.jld"

    # Save the array to JLD2 file
    @save jld_file_path num_strains_array
    
    return(nothing)
end

@time calc_num_strains()
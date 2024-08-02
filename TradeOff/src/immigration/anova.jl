using JLD2
using SimpleANOVA

function perform_anova()
    # Preallocate the variables I want to extract from the input
    num_immigrants = 0
    rl = 0
    ru = 0

    # Check that all arguments can be converted to integers
    try
        num_immigrants = parse(Int64, ARGS[1])
        rl = parse(Int64, ARGS[2])
        ru = parse(Int64, ARGS[3])
    catch e
        error("Need to provide an integer")
    end

    println("Compiled and input read in!")
    flush(stdout)

    frequencies = [10, 20, 40, 80, 160, 320, 640]

    # Initialize arrays to collect data
    community_EUE_array = []

    for i in frequencies
        # Construct the data directory and stats file paths
        data_dir = joinpath(pwd(), "Output", "niche_size$(rl)_$(ru)", "$(num_immigrants)immigrants", "$(i)events")
        stats_file = joinpath(data_dir, "RunStats$(i)events_$(num_immigrants)immigrants.jld")
        
        # Check it actually exists
        if !isfile(stats_file)
            error("missing stats file for $(i)events_$(num_immigrants)immigrants simulations")
        end

        # Load simulation data
        community_EUE = load(stats_file, "mean_community_EUE")

        # Collect data
        push!(community_EUE_array, community_EUE)
    end

    # Combine all data into a single vector and create a corresponding group vector
    combined_data = vcat(community_EUE_array...)
    group_labels = repeat(1:length(frequencies), inner=length(community_EUE_array[1]))

    # Perform the ANOVA test
    anova_result = anova(combined_data, group_labels)

    # Inspect the results
    println(anova_result)
end

@time perform_anova()

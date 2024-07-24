data_dir = joinpath(
    pwd(), "Output", "10immigrants", "10events")

    stats_file = joinpath(data_dir,  "AvRun$(rps)Data.jld")
    
    # Check it actually exists
    if ~isfile(stats_file)
        error("missing stats file for $(num_immigrations)events_$(num_immigrants)immigrants simulations")
    end

    # Extract time and surviving species data
    t_times = load(stats_file, "T")
    species_EUES = load(stats_file, "species_EUEs")
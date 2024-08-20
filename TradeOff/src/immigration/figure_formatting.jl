using Plots, JLD2

function figure1()
    frequencies = [10, 20, 40, 80, 160, 320, 640]
    # plot directory
    outdir1_5 = joinpath(pwd(), "Output", "niche_size1_5", "Immigration_plots")
    outdir20_25 = joinpath(pwd(), "Output", "niche_size20_25", "Immigration_plots")
    outdir1_25 = joinpath(pwd(), "Output", "niche_size1_25", "Immigration_plots")

    EUE_1_5 = load(joinpath(outdir1_5, "EUE_$(frequencies[1])to$(frequencies[end])_frequencies.png"))
    shannon_1_5 = load(joinpath(outdir1_5, "shannon_$(frequencies[1])to$(frequencies[end])_frequencies.png"))
    substrate_1_5 = load(joinpath(outdir1_5, "num_substrates_$(frequencies[1])to$(frequencies[end])_frequencies.png"))
    biomass_1_5 = load(joinpath(outdir1_5, "biomass_$(frequencies[1])to$(frequencies[end])_frequencies.png"))

    EUE_20_25 = load(joinpath(outdir20_25, "EUE_$(frequencies[1])to$(frequencies[end])_frequencies.png"))
    shannon_20_25 = load(joinpath(outdir20_25, "shannon_$(frequencies[1])to$(frequencies[end])_frequencies.png"))
    substrate_20_25 = load(joinpath(outdir20_25, "num_substrates_$(frequencies[1])to$(frequencies[end])_frequencies.png"))
    biomass_20_25 = load(joinpath(outdir20_25, "biomass_$(frequencies[1])to$(frequencies[end])_frequencies.png"))

    EUE_1_25 = load(joinpath(outdir1_25, "EUE_$(frequencies[1])to$(frequencies[end])_frequencies.png"))
    shannon_1_25 = load(joinpath(outdir1_25, "shannon_$(frequencies[1])to$(frequencies[end])_frequencies.png"))
    substrate_1_25 = load(joinpath(outdir1_25, "num_substrates_$(frequencies[1])to$(frequencies[end])_frequencies.png"))
    biomass_1_25 = load(joinpath(outdir1_25, "biomass_$(frequencies[1])to$(frequencies[end])_frequencies.png"))

    custom_layout = @layout [a b c d; e f g]

    figure= plot(
        EUE_1_5, EUE_20_25, EUE_1_25,
        shannon_1_5, shannon_20_25, shannon_1_25,
        substrate_1_5, substrate_20_25, substrate_1_25,
        biomass_1_5, biomass_20_25, biomass_1_25, 
        layout = (3, 4),
        size = (1800, 1600),
        margin = 10Plots.mm
    )

    # Define output directory and if necessary make it
    outdir = joinpath(pwd(), "Output", "figures")
    mkpath(outdir)

    savefig(figure, joinpath(outdir, "figure1.png"))
end

@time figure1()
using Plots
using TradeOff

function wrap_test()
    ωs, λlf, λhf, λlv, λhv = ω_test()
    # Define directory and if necessary make it
    plot_dir = joinpath(pwd(), "Output", "Fig8")
    mkpath(plot_dir)
    # Extract max ribosome fraction
    plot(ωs, λlf, label = "")
    savefig(joinpath(plot_dir, "GrowthLowFixed.png"))
    plot(ωs, λhf, label = "")
    savefig(joinpath(plot_dir, "GrowthHighFixed.png"))
    plot(ωs, λlv, label = "")
    savefig(joinpath(plot_dir, "GrowthLowVar.png"))
    plot(ωs, λhv, label = "")
    savefig(joinpath(plot_dir, "GrowthHighVar.png"))
    return (nothing)
end

@time wrap_test()

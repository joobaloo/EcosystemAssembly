using Plots, KernelDensity, JLD2, HypothesisTests
data = randn(10^4);

k = kde(data);
plot(k.x, k.density)

rl = 1
ru = 5
jld_file_path = "../Output/niche_size$(rl)_$(ru)/plateau_indices_array.jld2"
@load jld_file_path plateau_indices_array

println(plateau_indices_array[1])


# Example simplified data setup
EUE_data1_5 = [rand(1000) for _ in 1:8]
EUE_data20_25 = [rand(1000) for _ in 1:8]
EUE_data1_25 = [rand(1000) for _ in 1:8]
frequencies = 1:8

# Define color palette
custom_colour_pallete = [:lightblue, :darkblue, :red]

# Create a plot with the correct layout
distribution_plot = plot(layout = (length(frequencies), 1))

# Loop through each frequency and plot
for (i, freq) in enumerate(frequencies)
    hist_data1_5 = EUE_data1_5[i]
    hist_data20_25 = EUE_data20_25[i]
    hist_data1_25 = EUE_data1_25[i]

    # Concatenate the data for each frequency
    freq_data = hcat(hist_data1_5, hist_data20_25, hist_data1_25)

    # Print the distribution spread
    println(string(freq, " distribution spread is: ", maximum(freq_data) - minimum(freq_data), " from ",  minimum(freq_data), " to ", maximum(freq_data)))

    # Plot the histogram
    histogram!(
        distribution_plot, 
        freq_data, 
        bins=10, 
        fillalpha=0.5,
        label=["1-5", "20-25", "1-25"],
        color=custom_colour_pallete,
        title="Frequency $freq", 
        xlims=(0, 4),
        ylims=(0, 300),
        xlabel="EUE", 
        ylabel="Frequency",
        subplot=i
    )
end

# Display the plot
display(distribution_plot)

# Sample data
x = rand(100)  # Random data for x
y = 0.5*x + 0.1*rand(100)  # Data for y, with some noise

# Perform Pearson correlation test
result = CorrelationTest(x, y)

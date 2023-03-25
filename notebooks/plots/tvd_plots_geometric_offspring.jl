module TVDPlotsGeometricOffspring

using Pkg
Pkg.activate(".")

import TVDAlgorithm
import JSON
import Unmarshal
using Plots
using Plots.PlotMeasures

struct ConsolidatedTVDSimulationResults
    K_vals::Vector{Int}
    path_lengths::Vector{Int}
    num_trials::Int
    results::Matrix{TVDAlgorithm.TVDSimulationSummary}
end

file = open("out/data/tvd_geometric_offspring_z0_is_K.json", "r")
data = read(file, String)
close(file)

tvd_results_z0_is_K = Unmarshal.unmarshal(
    ConsolidatedTVDSimulationResults,
    JSON.parse(data)
)

file = open("out/data/tvd_geometric_offspring_z0_is_1.json", "r")
data = read(file, String)
close(file)

tvd_results_z0_is_1 = Unmarshal.unmarshal(
    ConsolidatedTVDSimulationResults,
    JSON.parse(data)
)

file = open("out/data/tvd_geometric_offspring_mean_only_z0_is_1.json", "r")
data = read(file, String)
close(file)

tvd_results_mean_only_z0_is_1 = Unmarshal.unmarshal(
    ConsolidatedTVDSimulationResults,
    JSON.parse(data)
)

path_length_plot = plot(
    1:50,
    i -> tvd_results_z0_is_1.results[i, 10].mean,
    title="Estimated TVD over different path lengths, for K=100 and z₀=1",
    xlabel="Path length",
    ylabel="Estimated TVD",
    linecolor=:purple,
    linealpha=0.9,
    linewidth=2,
    legend=false,
    titlefontsize=24,
    tickfontsize=16,
    guidefontsize=16,
    size=(1400, 1000),
    ylims=(0, 0.06),
    bottom_margin=12mm,
    left_margin=12mm,
    top_margin=3mm
)

savefig(path_length_plot, "out/plots/geometric_offspring/tvd_path_length.png")

comparative_path_length_plot = plot(
    1:50,
    [
        i -> tvd_results_z0_is_1.results[i, 10].mean
        i -> tvd_results_z0_is_K.results[i, 10].mean
    ],
    title="Estimated TVD over different path lengths, for K=100",
    xlabel="Path length",
    ylabel="Estimated TVD",
    linecolor=[:purple :thistle4],
    linealpha=0.9,
    linewidth=2,
    legend=:topleft,
    labels=["z₀ = 1" "z₀ = 100"],
    legendfontsize=16,
    titlefontsize=24,
    tickfontsize=16,
    guidefontsize=16,
    size=(1400, 1000),
    ylims=(0, 0.06),
    bottom_margin=12mm,
    left_margin=12mm,
    top_margin=3mm
)

savefig(comparative_path_length_plot, "out/plots/geometric_offspring/tvd_path_length_comparative.png")

mean_only_path_length_plot = plot(
    1:50,
    [
        i -> tvd_results_z0_is_1.results[i, 10].mean
        i -> tvd_results_mean_only_z0_is_1.results[i, 10].mean
    ],
    title="Estimated TVD over different path lengths, for K=100",
    xlabel="Path length",
    ylabel="Estimated TVD",
    linecolor=[:purple :thistle4],
    linealpha=0.9,
    linewidth=2,
    legend=:topleft,
    labels=["mean and variance matching" "mean matching only"],
    legendfontsize=16,
    titlefontsize=24,
    tickfontsize=16,
    guidefontsize=16,
    size=(1400, 1000),
    ylims=(0, 0.5),
    bottom_margin=12mm,
    left_margin=12mm,
    top_margin=3mm
)

savefig(mean_only_path_length_plot, "out/plots/geometric_offspring/tvd_path_length_mean_only.png")

K_plot = plot(
    i -> tvd_results_z0_is_1.K_vals[i],
    i -> tvd_results_z0_is_1.results[30, i].mean,
    1:20,
    title="Estimated TVD over different K, for path length 30 and z₀=1",
    xlabel="Carrying capacity",
    ylabel="Estimated TVD",
    linecolor=:purple,
    linealpha=0.9,
    linewidth=2,
    legend=false,
    titlefontsize=24,
    tickfontsize=16,
    guidefontsize=16,
    size=(1400, 1000),
    xticks=0:20:200,
    ylims=(0, 0.15),
    fillcolour=:thistle,
    fillalpha=0.4,
    bottom_margin=12mm,
    left_margin=12mm,
    top_margin=3mm
)

savefig(K_plot, "out/plots/geometric_offspring/tvd_K.png")

comparative_K_plot = plot(
    i -> tvd_results_z0_is_1.K_vals[i],
    [
        i -> tvd_results_z0_is_1.results[30, i].mean
        i -> tvd_results_z0_is_K.results[30, i].mean
    ],
    1:20,
    title="Estimated TVD over different K, for path length 30",
    xlabel="Carrying capacity",
    ylabel="Estimated TVD",
    linecolor=[:purple :thistle4],
    linealpha=0.9,
    linewidth=2,
    legend=:topright,
    labels=["z₀ = 1" "z₀ = 100"],
    legendfontsize=16,
    titlefontsize=24,
    tickfontsize=16,
    guidefontsize=16,
    size=(1400, 1000),
    ylims=(0, 0.15),
    bottom_margin=12mm,
    left_margin=12mm,
    top_margin=3mm
)

savefig(comparative_K_plot, "out/plots/geometric_offspring/tvd_K_comparative.png")

mean_only_K_plot = plot(
    i -> tvd_results_z0_is_1.K_vals[i],
    [
        i -> tvd_results_z0_is_1.results[30, i].mean
        i -> tvd_results_mean_only_z0_is_1.results[30, i].mean
    ],
    1:20,
    title="Estimated TVD over different K, for path length 30",
    xlabel="Carrying capacity",
    ylabel="Estimated TVD",
    linecolor=[:purple :thistle4],
    linealpha=0.9,
    linewidth=2,
    legend=:topright,
    labels=["mean and variance matching" "mean matching only"],
    legendfontsize=16,
    titlefontsize=24,
    tickfontsize=16,
    guidefontsize=16,
    size=(1400, 1000),
    ylims=(0, 0.5),
    bottom_margin=12mm,
    left_margin=12mm,
    top_margin=3mm
)

savefig(mean_only_K_plot, "out/plots/geometric_offspring/tvd_K_mean_only.png")

K_plot_multiple_path_lengths = plot(
    i -> tvd_results_z0_is_1.K_vals[i],
    [
        i -> tvd_results_z0_is_1.results[5, i].mean,
        i -> tvd_results_z0_is_1.results[10, i].mean,
        i -> tvd_results_z0_is_1.results[15, i].mean,
        i -> tvd_results_z0_is_1.results[20, i].mean,
        i -> tvd_results_z0_is_1.results[25, i].mean,
        i -> tvd_results_z0_is_1.results[30, i].mean,
        i -> tvd_results_z0_is_1.results[35, i].mean,
    ],
    1:20,
    title="Estimated TVD over different K, for z₀=1 and different path lengths",
    xlabel="Carrying capacity",
    ylabel="Estimated TVD",
    linealpha=0.9,
    linewidth=2,
    legend=:topright,
    labels=hcat("path length 5", "path length 10", "path length 15", "path length 20",
        "path length 25", "path length 30", "path length 35", "path length 40"),
    legendfontsize=16,
    titlefontsize=24,
    tickfontsize=16,
    guidefontsize=16,
    size=(1400, 1000),
    xticks=0:20:200,
    ylims=(0, 0.2),
    palette=palette([:purple, :lavender, :skyblue, :olivedrab], 8),
    fillalpha=0.4,
    bottom_margin=12mm,
    left_margin=12mm,
    top_margin=3mm
)

savefig(K_plot_multiple_path_lengths, "out/plots/geometric_offspring/tvd_K_and_path_length.png")

end

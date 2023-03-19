module TVDPathLengthPlot

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

file = open("out/data/tvd_results_z0_is_1.json", "r")
data = read(file, String)
close(file)

tvd_results_z0_is_1 = Unmarshal.unmarshal(
    ConsolidatedTVDSimulationResults,
    JSON.parse(data)
)

file = open("out/data/tvd_results_z0_is_K.json", "r")
data = read(file, String)
close(file)

tvd_results_z0_is_K = Unmarshal.unmarshal(
    ConsolidatedTVDSimulationResults,
    JSON.parse(data)
)

TVD1(i) = tvd_results_z0_is_K.results[30, i].mean
TVDK(i) = tvd_results_z0_is_1.results[30, i].mean

K_plot = plot(
    i -> tvd_results_z0_is_1.K_vals[i],
    TVDK,
    5:29,
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
    ylims=(0, 0.12),
    fillcolour=:thistle,
    fillalpha=0.4,
    bottom_margin=12mm,
    left_margin=12mm,
    top_margin=3mm
)

savefig(K_plot, "out/plots/tvd_K_plot.png")

comparative_K_plot = plot(
    i -> tvd_results_z0_is_1.K_vals[i],
    [
        TVD1,
        TVDK
    ],
    5:29,
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
    ylims=(0, 0.12),
    bottom_margin=12mm,
    left_margin=12mm,
    top_margin=3mm
)

savefig(comparative_K_plot, "out/plots/tvd_K_plot_comparative.png")

K_plot_multiple_path_lengths = plot(
    i -> tvd_results_z0_is_1.K_vals[i],
    [
        i -> tvd_results_z0_is_K.results[5, i].mean,
        i -> tvd_results_z0_is_K.results[10, i].mean,
        i -> tvd_results_z0_is_K.results[15, i].mean,
        i -> tvd_results_z0_is_K.results[20, i].mean,
        i -> tvd_results_z0_is_K.results[25, i].mean,
        i -> tvd_results_z0_is_K.results[30, i].mean,
        i -> tvd_results_z0_is_K.results[35, i].mean,
        i -> tvd_results_z0_is_K.results[40, i].mean,
    ],
    5:29,
    title="Estimated TVD over different K, for z₀=1 and different path lengths",
    xlabel="Carrying capacity",
    ylabel="Estimated TVD",
    linealpha=0.9,
    linewidth=2,
    legend=:topright,
    labels=["path length 5" "path length 10" "path length 15" "path length 20" "path length 25" "path length 30" "path length 35" "path length 40"],
    legendfontsize=16,
    titlefontsize=24,
    tickfontsize=16,
    guidefontsize=16,
    size=(1400, 1000),
    xticks=0:20:200,
    ylims=(0, 0.12),
    palette=palette([:purple, :lavender, :skyblue, :olivedrab], 8),
    fillalpha=0.4,
    bottom_margin=12mm,
    left_margin=12mm,
    top_margin=3mm
)

savefig(K_plot_multiple_path_lengths, "out/plots/tvd_K_and_path_length_plot.png")

end

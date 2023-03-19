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

TVD1(i) = tvd_results_z0_is_K.results[i, 19].mean
TVDK(i) = tvd_results_z0_is_1.results[i, 19].mean

path_length_plot = plot(
    1:50,
    TVD1,
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
    ylims=(0, 0.03),
    bottom_margin=12mm,
    left_margin=12mm,
    top_margin=3mm
)

savefig(path_length_plot, "out/plots/tvd_path_length_plot.png")

comparative_path_length_plot = plot(
    1:50,
    [
        TVD1,
        TVDK
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
    ylims=(0, 0.03),
    bottom_margin=12mm,
    left_margin=12mm,
    top_margin=3mm
)

savefig(comparative_path_length_plot, "out/plots/tvd_path_length_plot_comparative.png")

end

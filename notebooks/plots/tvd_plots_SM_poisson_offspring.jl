module TVDPlotsSMPoissonOffspring

using Pkg
Pkg.activate(".")

import TVDAlgorithm
import JSON
import Unmarshal
using Plots
using Plots.PlotMeasures

struct ConsolidatedTVDSimulationResults
    z₀_list::Vector{Int}
    path_lengths::Vector{Int}
    num_trials::Int
    results::Matrix{TVDAlgorithm.TVDSimulationSummary}
end

file = open("out/data/tvd_sm_poisson_offspring_over_z0.json", "r")
data = read(file, String)
close(file)

tvd_results = Unmarshal.unmarshal(
    ConsolidatedTVDSimulationResults,
    JSON.parse(data)
)

path_length_plot = plot(
    1:50,
    i -> tvd_results.results[i, 1].mean,
    title="Estimated TVD over different path lengths, for z₀ = 1",
    xlabel="Path length",
    ylabel="Estimated TVD",
    linecolor=:purple,
    linealpha=0.9,
    linewidth=3,
    legend=false,
    titlefontsize=48,
    tickfontsize=32,
    guidefontsize=42,
    size=(2800, 2000),
    ylims=(0, 0.2),
    bottom_margin=28mm,
    left_margin=28mm,
    top_margin=6mm
)

savefig(path_length_plot, "out/plots/poisson_offspring_SM/tvd_path_length.png")

comparative_path_length_plot = plot(
    1:50,
    [
        i -> tvd_results.results[i, 1].mean,
        i -> tvd_results.results[i, 2].mean,
        i -> tvd_results.results[i, 3].mean,
        i -> tvd_results.results[i, 4].mean,
        i -> tvd_results.results[i, 5].mean,
        i -> tvd_results.results[i, 10].mean,
        i -> tvd_results.results[i, 12].mean,
        i -> tvd_results.results[i, 13].mean,
        i -> tvd_results.results[i, 14].mean,
        i -> tvd_results.results[i, 15].mean
    ],
    title="Estimated TVD over different path lengths, for various z₀",
    xlabel="Path length",
    ylabel="Estimated TVD",
    linealpha=0.9,
    linewidth=3,
    legend=:bottomright,
    labels=hcat("z₀ = 1", "z₀ = 2", "z₀ = 3", "z₀ = 4", "z₀ = 5", "z₀ = 10",
        "z₀ = 20", "z₀ = 30", "z₀ = 40", "z₀ = 50"),
    legendfontsize=32,
    titlefontsize=48,
    tickfontsize=32,
    guidefontsize=42,
    size=(2800, 2000),
    xticks=0:20:200,
    ylims=(0, 0.2),
    palette=palette([:purple, :thistle, :skyblue, :olivedrab], 10),
    fillalpha=0.4,
    bottom_margin=28mm,
    left_margin=28mm,
    top_margin=6mm
)

savefig(comparative_path_length_plot, "out/plots/poisson_offspring_SM/tvd_path_length_comparative.png")

z₀_plot = plot(
    i -> tvd_results.z₀_list[i],
    i -> tvd_results.results[1, i].mean,
    1:30,
    title="Estimated TVD over different z₀, for path length 1",
    xlabel="z₀",
    ylabel="Estimated TVD",
    linecolor=:purple,
    linealpha=0.9,
    linewidth=3,
    legend=false,
    titlefontsize=48,
    tickfontsize=32,
    guidefontsize=42,
    size=(2800, 2000),
    ylims=(0, 0.1),
    bottom_margin=28mm,
    left_margin=28mm,
    top_margin=6mm
)

savefig(path_length_plot, "out/plots/poisson_offspring_SM/tvd_z0.png")

comparative_z₀_plot = plot(
    i -> tvd_results.z₀_list[i],
    [
        i -> tvd_results.results[1, i].mean,
        i -> tvd_results.results[5, i].mean,
        i -> tvd_results.results[10, i].mean,
        i -> tvd_results.results[15, i].mean,
        i -> tvd_results.results[20, i].mean,
        i -> tvd_results.results[25, i].mean,
        i -> tvd_results.results[30, i].mean,
        i -> tvd_results.results[35, i].mean,
        i -> tvd_results.results[40, i].mean
    ],
    1:30,
    title="Estimated TVD over different z₀, for various path lengths",
    xlabel="z₀",
    ylabel="Estimated TVD",
    linealpha=0.9,
    linewidth=3,
    legend=:topright,
    labels=hcat("path length 1", "path length 5", "path length 10", "path length 15", "path length 20",
        "path length 25", "path length 30", "path length 35", "path length 40"),
    legendfontsize=32,
    titlefontsize=48,
    tickfontsize=32,
    guidefontsize=42,
    size=(2800, 2000),
    xticks=0:20:200,
    ylims=(0, 0.2),
    palette=palette([:purple, :thistle, :skyblue, :olivedrab], 9),
    fillalpha=0.4,
    bottom_margin=28mm,
    left_margin=28mm,
    top_margin=6mm
)

savefig(comparative_z₀_plot, "out/plots/poisson_offspring_SM/tvd_z0_comparative.png")

end
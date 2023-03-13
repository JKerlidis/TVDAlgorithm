module TVDSurfacePlot

using Pkg
Pkg.activate(".")

import TVDAlgorithm
import JSON
import Unmarshal
using Plots
using Plots.PlotMeasures

file = open("data/tvd_data.json", "r")
data = read(file, String)
close(file)

results = Unmarshal.unmarshal(
    Array{Array{TVDAlgorithm.TVDSimulationSummary}},
    JSON.parse(data)
)
results = hcat(results...)

path_length_plot = plot(
    path_lengths,
    i -> results[i, 10].mean,
    title="Estimated TVD over different path lengths, for K=100",
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
    ylims=(0, 0.1),
    ribbon=[results[i, 10].standard_deviation for i ∈ path_lengths],
    fillcolour=:thistle,
    fillalpha=0.4,
    bottom_margin=12mm,
    left_margin=12mm,
    top_margin=3mm
)

savefig(path_length_plot, "data/tvd_path_length_plot.png")

K_plot = plot(
    1:20,
    i -> results[10, i].mean,
    title="Estimated TVD over different K, for path length 10",
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
    xformatter=i -> floor(10 * i),
    ylims=(0, 0.1),
    ribbon=[results[10, i].standard_deviation for i ∈ 1:20],
    fillcolour=:thistle,
    fillalpha=0.4,
    bottom_margin=12mm,
    left_margin=12mm,
    top_margin=3mm
)

savefig(K_plot, "data/tvd_K_plot.png")

end
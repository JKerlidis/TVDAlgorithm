module TVDPlotsPoissonOffspring

using Pkg
Pkg.activate(".")

import TVDAlgorithm
import JSON
import Unmarshal
using Plots
using Plots.PlotMeasures

struct ZeroInflationTVDResults
    K::Int
    M::Int
    num_trials::Int
    path_length::Int
    z0_values::Vector{Int}
    lambda_values::Vector{Float64}
    results::Matrix{Float64}
end

struct ConsolidatedTVDSimulationResults
    K_vals::Vector{Int}
    path_lengths::Vector{Int}
    num_trials::Int
    results::Matrix{TVDAlgorithm.TVDSimulationSummary}
end

file = open("out/data/tvd_poisson_offspring_zero_inflation.json", "r")
data = read(file, String)
close(file)

zero_inflation_tvd_results = Unmarshal.unmarshal(
    ZeroInflationTVDResults,
    JSON.parse(data)
)

file = open("out/data/tvd_poisson_offspring_z0_is_K.json", "r")
data = read(file, String)
close(file)

tvd_results_z0_is_K = Unmarshal.unmarshal(
    ConsolidatedTVDSimulationResults,
    JSON.parse(data)
)

file = open("out/data/tvd_poisson_offspring_z0_is_1.json", "r")
data = read(file, String)
close(file)

tvd_results_z0_is_1 = Unmarshal.unmarshal(
    ConsolidatedTVDSimulationResults,
    JSON.parse(data)
)

file = open("out/data/tvd_poisson_offspring_mean_only_z0_is_1.json", "r")
data = read(file, String)
close(file)

tvd_results_mean_only_z0_is_1 = Unmarshal.unmarshal(
    ConsolidatedTVDSimulationResults,
    JSON.parse(data)
)


zero_inflation_plot = plot(
    1:100,
    i -> zero_inflation_tvd_results.results[6, i],
    title="Estimated one-step TVD over z₀, for K=100, M=2, and λ=2.5",
    xlabel="Initial population size (z₀)",
    ylabel="Estimated TVD",
    linecolor=:purple,
    linealpha=0.9,
    linewidth=2,
    legend=false,
    titlefontsize=24,
    tickfontsize=16,
    guidefontsize=16,
    size=(1400, 1000),
    ylims=(0, 0.05),
    bottom_margin=12mm,
    left_margin=12mm,
    top_margin=3mm
)

savefig(zero_inflation_plot, "out/plots/poisson_offspring/zero_inflation.png")

zero_inflation_comparative_plot = plot(
    1:100,
    [
        i -> zero_inflation_tvd_results.results[1, i],
        i -> zero_inflation_tvd_results.results[2, i],
        i -> zero_inflation_tvd_results.results[3, i],
        i -> zero_inflation_tvd_results.results[4, i],
        i -> zero_inflation_tvd_results.results[5, i],
        i -> zero_inflation_tvd_results.results[6, i],
        i -> zero_inflation_tvd_results.results[7, i],
        i -> zero_inflation_tvd_results.results[8, i],
        i -> zero_inflation_tvd_results.results[9, i],
        i -> zero_inflation_tvd_results.results[10, i],
        i -> zero_inflation_tvd_results.results[11, i],
    ],
    title="Estimated one-step TVD over z₀, for K=100, M=2, and various λ",
    xlabel="Initial population size (z₀)",
    ylabel="Estimated TVD",
    linealpha=0.9,
    linewidth=2,
    legend=:topright,
    labels=hcat("λ=2.0", "λ=2.1", "λ=2.2", "λ=2.3", "λ=2.4", "λ=2.5",
        "λ=2.6", "λ=2.7", "λ=2.8", "λ=2.9", "λ=3.0"),
    legendfontsize=16,
    titlefontsize=24,
    tickfontsize=16,
    guidefontsize=16,
    size=(1400, 1000),
    xticks=0:20:200,
    ylims=(0, 0.08),
    palette=palette([:purple, :lavender, :skyblue, :olivedrab], 11),
    fillalpha=0.4,
    bottom_margin=12mm,
    left_margin=12mm,
    top_margin=3mm
)

savefig(zero_inflation_comparative_plot, "out/plots/poisson_offspring/zero_inflation_comparative.png")

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
    ylims=(0, 0.14),
    bottom_margin=12mm,
    left_margin=12mm,
    top_margin=3mm
)

savefig(path_length_plot, "out/plots/poisson_offspring/tvd_path_length.png")

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
    ylims=(0, 0.14),
    bottom_margin=12mm,
    left_margin=12mm,
    top_margin=3mm
)

savefig(comparative_path_length_plot, "out/plots/poisson_offspring/tvd_path_length_comparative.png")

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
    ylims=(0, 0.7),
    bottom_margin=12mm,
    left_margin=12mm,
    top_margin=3mm
)

savefig(mean_only_path_length_plot, "out/plots/poisson_offspring/tvd_path_length_mean_only.png")

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
    ylims=(0, 0.35),
    fillcolour=:thistle,
    fillalpha=0.4,
    bottom_margin=12mm,
    left_margin=12mm,
    top_margin=3mm
)

savefig(K_plot, "out/plots/poisson_offspring/tvd_K.png")

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
    labels=["z₀ = 1" "z₀ = K"],
    legendfontsize=16,
    titlefontsize=24,
    tickfontsize=16,
    guidefontsize=16,
    size=(1400, 1000),
    ylims=(0, 0.35),
    bottom_margin=12mm,
    left_margin=12mm,
    top_margin=3mm
)

savefig(comparative_K_plot, "out/plots/poisson_offspring/tvd_K_comparative.png")

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
    ylims=(0, 0.8),
    bottom_margin=12mm,
    left_margin=12mm,
    top_margin=3mm
)

savefig(mean_only_K_plot, "out/plots/poisson_offspring/tvd_K_mean_only.png")

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
        i -> tvd_results_z0_is_1.results[40, i].mean,
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
    ylims=(0, 0.4),
    palette=palette([:purple, :lavender, :skyblue, :olivedrab], 8),
    fillalpha=0.4,
    bottom_margin=12mm,
    left_margin=12mm,
    top_margin=3mm
)

savefig(K_plot_multiple_path_lengths, "out/plots/poisson_offspring/tvd_K_and_path_length.png")

end

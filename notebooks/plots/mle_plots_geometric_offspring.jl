module MLEPlotsGeometricOffspring

using Pkg
Pkg.activate(".")

import TVDAlgorithm
import JSON
import Unmarshal
using Plots
using Plots.PlotMeasures

struct MLEResults
    K_vals::Vector{Int}
    path_length::Int
    num_trials::Int
    true_model::String
    alternate_model::String
    results::Vector{Int}
end

file = open("out/data/mle_geometric_offspring_psdbp_true.json")
data = read(file, String)
close(file)

mle_results_psdbp_true = Unmarshal.unmarshal(
    MLEResults,
    JSON.parse(data)
)

file = open("out/data/mle_geometric_offspring_cbp_true.json")
data = read(file, String)
close(file)

mle_results_cbp_true = Unmarshal.unmarshal(
    MLEResults,
    JSON.parse(data)
)

psdbp_p(i) = mle_results_psdbp_true.results[i] / mle_results_psdbp_true.num_trials
cbp_p(i) = mle_results_cbp_true.results[i] / mle_results_cbp_true.num_trials
psdbp_error = [
    1.96 * √(psdbp_p(i) * (1 - psdbp_p(i)) / mle_results_psdbp_true.num_trials)
    for i ∈ 1:length(mle_results_psdbp_true.K_vals)
]
cbp_error = [
    1.96 * √(cbp_p(i) * (1 - cbp_p(i)) / mle_results_cbp_true.num_trials)
    for i ∈ 1:length(mle_results_cbp_true.K_vals)
]

mle_K_plot = plot(
    i -> mle_results_psdbp_true.K_vals[i],
    psdbp_p,
    1:length(mle_results_psdbp_true.K_vals),
    ribbon=psdbp_error,
    title="Accuracy in selecting the generating model of a path, for various K",
    xlabel="Carrying Capacity",
    ylabel="Accuracy",
    linecolor=:purple,
    linealpha=0.9,
    linewidth=2,
    fillcolor=:purple,
    fillalpha=0.2,
    legend=:topright,
    label="PSDBP true model",
    legendfontsize=16,
    titlefontsize=24,
    tickfontsize=16,
    guidefontsize=16,
    size=(1400, 1000),
    ylims=(0, 1),
    bottom_margin=12mm,
    left_margin=12mm,
    top_margin=3mm
)

mle_K_plot = plot!(
    i -> mle_results_cbp_true.K_vals[i],
    cbp_p,
    1:length(mle_results_cbp_true.K_vals),
    ribbon=cbp_error,
    linecolor=:thistle4,
    linealpha=0.9,
    linewidth=2,
    fillcolor=:thistle4,
    fillalpha=0.2,
    label="CBP true model"
)

mle_K_plot = hline!(
    [0.5],
    linecolor=:grey,
    linealpha=0.4,
    label=false
)

savefig(mle_K_plot, "out/plots/geometric_offspring/MLE_K.png")

end

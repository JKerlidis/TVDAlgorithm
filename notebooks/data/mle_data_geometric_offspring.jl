module MLEDataGeometricOffspring

using Pkg
Pkg.activate(".")

import TVDAlgorithm
import Random
import ProgressMeter
import JSON

struct MLEModel
    model::Type
    K::Integer
end

rng = Random.Xoshiro(613)
K_vals = [i for i ∈ 100:100:150]
path_length = 30
num_trials = 10

# For various values of K, simulate from a true distribution, and then find the
# MLE on both the true distribution and an alternate distribution, in
# particular whether the likelihood is better maximised by the true model than
# the alternate
function run_simulation(
    rng::Random.AbstractRNG,
    K_vals::Vector{<:Integer},
    path_length::Integer,
    num_trials::Integer,
    true_distribution::TVDAlgorithm.BranchingProcess,
    alternate_distribution::TVDAlgorithm.BranchingProcess,
    z₀::Integer,
    output_file::String
)
    num_correctly_identified = zeros(Int, length(K_vals))
    selected_models = Matrix{MLEModel}(undef, length(K_vals), num_trials)

    for (i, K) ∈ enumerate(K_vals)

        num_correct = 0
        ProgressMeter.@showprogress string("Running simulation for K=", K) for j ∈ 1:num_trials
            selected_model = TVDAlgorithm.select_model_in_K(
                rng,
                TVDAlgorithm.transition_probabilities(
                    TVDAlgorithm.substitute_K(true_distribution, K, true)
                ),
                true_distribution,
                alternate_distribution,
                path_length,
                z₀
            )
            selected_models[i, j] = MLEModel(
                typeof(selected_model),
                selected_model.K
            )
            if typeof(selected_model) == typeof(true_distribution)
                num_correct += 1
            end
        end
        num_correctly_identified[i] = num_correct
    end

    file = open(output_file, "w")
    write(
        file,
        JSON.json(Dict(
            "K_vals" => K_vals,
            "path_length" => path_length,
            "num_trials" => num_trials,
            "true_model" => typeof(true_distribution),
            "alternate_model" => typeof(alternate_distribution),
            "data" => selected_models,
            "results" => num_correctly_identified
        ))
    )
    close(file)

    return num_correctly_identified
end

println("Simulating paths from the PSDBP:")
psdbp_mle_num_correct = run_simulation(
    rng,
    K_vals,
    path_length,
    num_trials,
    TVDAlgorithm.PSDBPMatchingKGeometricOffspring(1, 2),
    TVDAlgorithm.CBPKGeometricOffspring(1, 2),
    1,
    "out/data/mle_geometric_offspring_psdbp_true.json"
)

println("Simulating paths from the CBP:")
cbp_mle_num_correct = run_simulation(
    rng,
    K_vals,
    path_length,
    num_trials,
    TVDAlgorithm.CBPKGeometricOffspring(1, 2),
    TVDAlgorithm.PSDBPMatchingKGeometricOffspring(1, 2),
    1,
    "out/data/mle_geometric_offspring_cbp_true.json"
)

println("PSDBP results: ", psdbp_mle_num_correct)
println("CBP results: ", cbp_mle_num_correct)

end

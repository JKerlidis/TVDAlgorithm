module Simulations

import TVDAlgorithm
import Random
import ProgressMeter
import JSON

# Simulate the TVD between the CBP and PSDBP over different carrying capacities
# (K) and different path lengths
function tvd_simulation(
    rng::Random.AbstractRNG,
    K_vals::Vector{<:Integer},
    path_length::Integer,
    num_trials::Integer,
    cbp_distribution::TVDAlgorithm.BranchingProcess,
    psdbp_distribution::TVDAlgorithm.BranchingProcess,
    first_generation_size::Union{Integer,String},
    output_file::String
)
    tvd_results = []

    for K ∈ K_vals
        println("Simulating TVD for K = ", K)

        z₀ = 1
        if typeof(first_generation_size) <: Integer
            z₀ = first_generation_size
        elseif first_generation_size == "K"
            z₀ = K
        end

        cbp_transition_probs = TVDAlgorithm.transition_probabilities(
            TVDAlgorithm.substitute_K(cbp_distribution, K, true)
        )
        psdbp_transition_probs = TVDAlgorithm.transition_probabilities(
            TVDAlgorithm.substitute_K(psdbp_distribution, K, true)
        )

        simulation_output = TVDAlgorithm.approximate_tvd_extended_output(
            rng,
            num_trials,
            z₀,
            path_length,
            psdbp_transition_probs,
            cbp_transition_probs
        )

        tvd_results = tvd_results == [] ? simulation_output : [tvd_results simulation_output]
    end

    file = open(output_file, "w")
    write(
        file,
        JSON.json(Dict(
            "K_vals" => K_vals,
            "path_lengths" => [i for i ∈ 1:path_length],
            "num_trials" => num_trials,
            "results" => tvd_results
        ))
    )
    close(file)
end

struct MLEModel{T<:Real}
    model::Type
    K::T
end

# For various values of K, simulate from a true distribution, and then find the
# MLE on both the true distribution and an alternate distribution, in
# particular whether the likelihood is better maximised by the true model than
# the alternate
function mle_simulation(
    rng::Random.AbstractRNG,
    K_vals::Vector{<:Real},
    path_length::Integer,
    num_trials::Integer,
    true_distribution::TVDAlgorithm.BranchingProcess,
    alternate_distribution::TVDAlgorithm.BranchingProcess,
    z₀::Integer,
    output_file::String
)

    num_correctly_identified = zeros(Int, length(K_vals))
    selected_models = Matrix{MLEModel{eltype(K_vals)}}(undef, length(K_vals), num_trials)

    for (i, K) ∈ enumerate(K_vals)

        num_correct = 0
        transition_probabilities = TVDAlgorithm.transition_probabilities(
            TVDAlgorithm.substitute_K(true_distribution, K, true)
        )

        ProgressMeter.@showprogress string("Running simulation for K=", K) for j ∈ 1:num_trials
            selected_model = TVDAlgorithm.select_model_in_K(
                rng,
                transition_probabilities,
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

end

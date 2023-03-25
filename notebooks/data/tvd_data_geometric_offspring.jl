module TVDDataGeometricOffspring

using Pkg
Pkg.activate(".")

import TVDAlgorithm
import Random
import JSON

rng = Random.Xoshiro(415)
K_vals = [i for i ∈ 10:10:200]
path_length = 50
num_trials = 1000000

# Simulate the TVD between the CBP and PSDBP over different carrying capacities
# (K) and different path lengths
function run_simulation(
    rng::Random.AbstractRNG,
    K_vals::Vector{<:Integer},
    path_length::Integer,
    num_trials::Integer,
    first_generation_size::Union{Integer,String},
    mean_only::Bool,
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
            TVDAlgorithm.CBPKGeometricOffspring(K)
        )
        psdbp_transition_probs = TVDAlgorithm.transition_probabilities(
            TVDAlgorithm.PSDBPMatchingKGeometricOffspring(K, mean_only)
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

    return tvd_results
end

run_simulation(
    rng, K_vals, path_length, num_trials, "K", false, "out/data/tvd_geometric_offspring_z0_is_K.json"
)

run_simulation(
    rng, K_vals, path_length, num_trials, 1, false, "out/data/tvd_geometric_offspring_z0_is_1.json"
)

run_simulation(
    rng, K_vals, path_length, num_trials, 1, true, "out/data/tvd_geometric_offspring_mean_only_z0_is_1.json"
)

end

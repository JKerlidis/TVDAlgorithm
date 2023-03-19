module TVDDataGeneration

using Pkg
Pkg.activate(".")

import TVDAlgorithm
import Random
import JSON


rng = Random.Xoshiro(514)
K_vals = [[i for i ∈ 1:9]; [i for i ∈ 10:10:200]]
path_length = 50
num_trials = 1000000

# Simulate the TVD between the CBP and PSDBP over different carrying capacities
# (K), different path lengths, and for z₀=1 and z₀=K
function run_simulation(
    K_vals::Vector{<:Integer},
    path_length::Integer,
    num_trials::Integer
)
    tvd_results_z₀_is_1 = []
    tvd_results_z₀_is_K = []

    for K ∈ K_vals
        println("Simulating TVD for K = ", K)

        cbp_transition_probs = TVDAlgorithm.transition_probabilities(
            TVDAlgorithm.CBPCarryingCapacityBinomial(K)
        )
        psdbp_transition_probs = TVDAlgorithm.transition_probabilities(
            TVDAlgorithm.PSDBPCarryingCapacityNegativeBinomial(K)
        )

        simulation_output_1 = TVDAlgorithm.approximate_tvd_extended_output(
            rng,
            num_trials,
            1,
            path_length,
            psdbp_transition_probs,
            cbp_transition_probs
        )

        if tvd_results_z₀_is_1 == []
            tvd_results_z₀_is_1 = simulation_output_1
        else
            tvd_results_z₀_is_1 = [tvd_results_z₀_is_1 simulation_output_1]
        end

        simulation_output_K = TVDAlgorithm.approximate_tvd_extended_output(
            rng,
            num_trials,
            K,
            path_length,
            psdbp_transition_probs,
            cbp_transition_probs
        )

        if tvd_results_z₀_is_K == []
            tvd_results_z₀_is_K = simulation_output_K
        else
            tvd_results_z₀_is_K = [tvd_results_z₀_is_K simulation_output_K]
        end
    end

    return tvd_results_z₀_is_1, tvd_results_z₀_is_K
end

tvd_results_z₀_is_1, tvd_results_z₀_is_K = run_simulation(
    K_vals, path_length, num_trials
)

# Save the results to file
file = open("out/data/tvd_results_z0_is_1.json", "w")
write(
    file,
    JSON.json(Dict(
        "K_vals" => K_vals,
        "path_lengths" => [i for i ∈ 1:path_length],
        "num_trials" => num_trials,
        "results" => tvd_results_z₀_is_1
    ))
)
close(file)

file = open("out/data/tvd_results_z0_is_K.json", "w")
write(
    file,
    JSON.json(Dict(
        "K_vals" => K_vals,
        "path_lengths" => [i for i ∈ 1:path_length],
        "num_trials" => num_trials,
        "results" => tvd_results_z₀_is_K
    ))
)
close(file)

end

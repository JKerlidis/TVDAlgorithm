module TVDDataGeneration

using Pkg
Pkg.activate(".")

import TVDAlgorithm
import Random
import JSON


rng = Random.Xoshiro(514)
K_vals = [i for i ∈ 10:10:200]
path_length = 50
num_trials = 1000000

results = []
for K ∈ K_vals
    println("Simulating K = ", K)
    cbp = TVDAlgorithm.CBPCarryingCapacityBinomial(K)
    psdbp = TVDAlgorithm.PSDBPCarryingCapacityNegativeBinomial(K)

    simulation_output = TVDAlgorithm.approximate_tvd_extended_output(
        rng,
        num_trials,
        K,
        path_length,
        transition_probabilities(psdbp),
        transition_probabilities(cbp)
    )

    if results == []
        results = simulation_output
    else
        results = [results simulation_output]
    end
end

file = open("data/tvd_data.json", "w")
write(file, JSON.json(results))
close(file)

end
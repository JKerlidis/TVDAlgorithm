module MLEIntGeometricOffspring

include("./simulations.jl")

import .Simulations
import TVDAlgorithm
import Random

K_vals = [i for i ∈ 10:10:200]
path_length = 30
num_trials = 5000

println("Simulating paths from the PSDBP:")
psdbp_mle_num_correct = Simulations.mle_simulation(
    Random.Xoshiro(613),
    K_vals,
    path_length,
    num_trials,
    TVDAlgorithm.PSDBPMatchingKGeometricOffspring(1, 2),
    TVDAlgorithm.CBPKGeometricOffspring(1, 2),
    1,
    "out/data/mle_int_geometric_offspring_psdbp_true.json"
)

println("Simulating paths from the CBP:")
cbp_mle_num_correct = Simulations.mle_simulation(
    Random.Xoshiro(614),
    K_vals,
    path_length,
    num_trials,
    TVDAlgorithm.CBPKGeometricOffspring(1, 2),
    TVDAlgorithm.PSDBPMatchingKGeometricOffspring(1, 2),
    1,
    "out/data/mle_int_geometric_offspring_cbp_true.json"
)

println("PSDBP results: ", psdbp_mle_num_correct)
println("CBP results: ", cbp_mle_num_correct)

end

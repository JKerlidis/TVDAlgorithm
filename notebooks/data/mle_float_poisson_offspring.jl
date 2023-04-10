module MLEFloatPoissonOffspring

include("./simulations.jl")

import .Simulations
import TVDAlgorithm
import Random

K_vals = [i for i ∈ 10.0:10:200.0]
path_length = 30
num_trials = 5000

println("Simulating paths from the PSDBP:")
psdbp_mle_num_correct = Simulations.mle_simulation(
    Random.Xoshiro(7539),
    K_vals,
    path_length,
    num_trials,
    TVDAlgorithm.PSDBPMatchingKPoissonOffspring(1.0, 2),
    TVDAlgorithm.CBPKPoissonOffspring(1.0, 2),
    1,
    "out/data/mle_float_poisson_offspring_psdbp_true.json"
)

println("Simulating paths from the CBP:")
cbp_mle_num_correct = Simulations.mle_simulation(
    Random.Xoshiro(7540),
    K_vals,
    path_length,
    num_trials,
    TVDAlgorithm.CBPKPoissonOffspring(1.0, 2),
    TVDAlgorithm.PSDBPMatchingKPoissonOffspring(1.0, 2),
    1,
    "out/data/mle_float_poisson_offspring_cbp_true.json"
)

println("PSDBP results: ", psdbp_mle_num_correct)
println("CBP results: ", cbp_mle_num_correct)

end

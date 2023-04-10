module TVDDataBinomialOffsring

include("./simulations.jl")

import .Simulations
import TVDAlgorithm
import Random

rng = Random.Xoshiro(514)
K_vals = [[i for i ∈ 5:9]; [i for i ∈ 10:10:200]]
path_length = 50
num_trials = 1000000

# TVD simulation with z₀ = 1
Simulations.tvd_simulation(
    rng,
    K_vals,
    path_length,
    num_trials,
    TVDAlgorithm.CBPKBinomialOffspring(1),
    TVDAlgorithm.PSDBPMatchingKBinomialOffspring(1),
    1,
    "out/data/tvd_binomial_offspring_z0_is_1.json"
)

# TVD simulation with z₀ = K
Simulations.tvd_simulation(
    rng,
    K_vals,
    path_length,
    num_trials,
    TVDAlgorithm.CBPKBinomialOffspring(1),
    TVDAlgorithm.PSDBPMatchingKBinomialOffspring(1),
    "K",
    "out/data/tvd_binomial_offspring_z0_is_K.json"
)

end

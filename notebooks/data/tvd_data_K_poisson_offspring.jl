module TVDDataKPoissonOffspring

include("./simulations.jl")

import .Simulations
import TVDAlgorithm
import Random

# Test the impact of changing z₀ and λ on the TVD
Simulations.zero_inflation_tvd_K_simulation(
    Random.Xoshiro(7465),
    100,
    2,
    1,
    1000000,
    [i for i ∈ 1:100],
    [i for i ∈ 2.0:0.1:3.0],
    "out/data/tvd_poisson_offspring_zero_inflation.json"
)

rng = Random.Xoshiro(415)
K_vals = [i for i ∈ 10:10:200]
path_length = 50
num_trials = 1000000

# TVD simulation with z₀ = K
Simulations.tvd_K_simulation(
    rng,
    K_vals,
    path_length,
    num_trials,
    TVDAlgorithm.CBPKPoissonOffspring(100, 2, 3.0),
    TVDAlgorithm.PSDBPMatchingKPoissonOffspring(100, 2, 3.0),
    "K",
    "out/data/tvd_poisson_offspring_z0_is_K.json"
)

# TVD simulation with z₀ = 1
Simulations.tvd_K_simulation(
    rng,
    K_vals,
    path_length,
    num_trials,
    TVDAlgorithm.CBPKPoissonOffspring(100, 2, 3.0),
    TVDAlgorithm.PSDBPMatchingKPoissonOffspring(100, 2, 3.0),
    1,
    "out/data/tvd_poisson_offspring_z0_is_1.json"
)

# TVD simulation with z₀ = 1 and processes matching in mean only
Simulations.tvd_K_simulation(
    rng,
    K_vals,
    path_length,
    num_trials,
    TVDAlgorithm.CBPKPoissonOffspring(100, 2, 3.0),
    TVDAlgorithm.PSDBPMatchingKPoissonOffspring(100, 2, 3.0, true),
    1,
    "out/data/tvd_poisson_offspring_mean_only_z0_is_1.json"
)

end

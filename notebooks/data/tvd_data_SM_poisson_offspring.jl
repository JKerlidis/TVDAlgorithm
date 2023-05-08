module TVDDataSMPoissonOffspring

include("./simulations.jl")

import .Simulations
import TVDAlgorithm
import Random

rng = Random.Xoshiro(292)
z₀_list = [[i for i ∈ 1:10]; [i for i ∈ 10:10:200]]
path_length = 50
num_trials = 1000000

# TVD simulation over z₀
Simulations.tvd_z₀_simulation(
    rng,
    z₀_list,
    path_length,
    num_trials,
    TVDAlgorithm.CBPSMPoissonOffspring(2, 2.0, 1000),
    TVDAlgorithm.PSDBPMatchingSMPoissonOffspring(2, 2.0, 1000),
    "out/data/tvd_sm_poisson_offspring_over_z0.json"
)

end

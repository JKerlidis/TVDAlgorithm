module TVDAlgorithm

import Distributions
import Random

export
    # Helper functions
    index_position,
    position_index,
    kahan_sum,

    # Exact TVD algorithm
    k_step_transition_probabilities,
    exact_tvd,

    # Approximate TVD algorithm
    CensoredObservation,
    simulate_coupling_probability,
    approximate_tvd

include("helpers.jl")
include("exact_tvd_algorithm.jl")
include("approximate_tvd_algorithm.jl")
include("processes.jl")

end

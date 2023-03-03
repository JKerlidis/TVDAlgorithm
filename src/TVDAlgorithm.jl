module TVDAlgorithm

using Random

export

    # Exact TVD algorithm
    index_position,
    position_index,
    k_step_transition_probabilities,
    kahan_sum,
    exact_tvd,

    # Approximate TVD algorithm
    CensoredObservation,
    simulate_coupling_probability,
    approximate_tvd

include("exact_tvd_algorithm.jl")
include("approximate_tvd_algorithm.jl")

end

module TVDAlgorithm

using Random

export
    # Compensated summation
    kahan_sum,

    # Exact TVD algorithm
    index_position,
    position_index,
    k_step_transition_probabilities,
    exact_tvd,

    # Approximate TVD algorithm
    simulate_coupling_probability

include("compensated_summation.jl")
include("exact_tvd_algorithm.jl")
include("approximate_tvd_algorithm.jl")

end
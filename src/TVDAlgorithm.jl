module TVDAlgorithm

using Distributions

export
    # Exact TVD algorithm
    index_position,
    position_index,
    kahan_sum,
    k_step_transition_probabilities,
    exact_tvd_algorithm

include("exact_tvd_algorithm.jl")

end

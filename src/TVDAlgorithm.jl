module TVDAlgorithm

import Distributions
import Random
import Optim

export
    # Helper functions
    index_position,
    position_index,
    CensoredObservation,

    # Branching processes
    BranchingProcess,
    CBPKBinomialOffspring,
    CBPKGeometricOffspring,
    PSDBPMatchingKBinomialOffspring,
    PSDBPMatchingKGeometricOffspring,
    substitute_K,
    transition_probabilities,
    log_likelihood,

    # Exact TVD algorithm
    k_step_transition_probabilities,
    exact_tvd,

    # Approximate TVD algorithm
    TVDSimulationOutput,
    TVDSimulationSummary,
    summarise,
    sample_path,
    simulate_coupling_probability,
    approximate_tvd,
    approximate_tvd_extended_output,

    # Model selection
    maximise_likelihood_in_K,
    select_model_in_K

include("helpers.jl")
include("processes.jl")
include("exact_tvd_algorithm.jl")
include("approximate_tvd_algorithm.jl")
include("model_selection.jl")

end

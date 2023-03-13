# Given a square matrix of transition probabilities for a homogeneous Markov
# chain, produce an (k+1)-tensor listing the probabilites of each path of
# length k (i.e. the k-step transition probabilities).
function k_step_transition_probabilities(
    k::Integer,
    P::Matrix{Float64}
)::Array{Float64}

    size(P, 1) ≠ size(P, 2) && throw(DimensionMismatch("P should be a square matrix"))
    k < 1 && throw(DomainError(k, "number of steps should be positive"))

    n = size(P, 1)
    d = k + 1

    k_step_probabilities = Array{Float64}(undef, ntuple(i -> size(P, 1), d)) # empty (k+1)-tensor

    for path_index ∈ eachindex(k_step_probabilities)
        path_position = index_position(path_index, n, d)

        path_probability = 1
        for t ∈ 1:k  # index of the t'th transition
            path_probability *= P[path_position[t], path_position[t+1]]
        end

        k_step_probabilities[path_index] = path_probability
    end

    k_step_probabilities
end

# Given two homogeneous Markov processes with the same finite state space and
# transition probability matrices P and Q respectively, calculate the total
# variation distance between them for paths of length k, for a given initial
# state z₀. The path space quickly becomes large as k increases, which can
# cause memory issues with large path lengths
function exact_tvd(
    z₀::Integer,
    k::Integer,
    P::Matrix{Float64},
    Q::Matrix{Float64}
)::Float64

    size(P) ≠ size(Q) && throw(DimensionMismatch("P and Q do not have the same state space"))
    size(P, 1) ≠ size(P, 2) && throw(DimensionMismatch("P and Q are not square matrices"))
    (z₀ < 0 || z₀ ≥ size(P, 1)) && throw(DomainError("z₀ should be in the range [0, √|P| - 1]"))

    n = size(P, 1)

    min_k_step_probs = min.(
        k_step_transition_probabilities(k, P),
        k_step_transition_probabilities(k, Q)
    )

    z₀_path_min_probs = min_k_step_probs[1:n:n^(k+1)]
    tvd = 1 - kahan_sum(z₀_path_min_probs)

    tvd
end

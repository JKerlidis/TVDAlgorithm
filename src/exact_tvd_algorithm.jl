# Convert the ith index of a d-tensor with each dimension containing n elements
# to a tuple representing its position in the tensor. If n==3 and d==3, then
# an index of 11 will return (2, 1, 2), because 11 is in the second row, first
# column, and at depth one in the 3-tensor.
function index_position(index::Integer, n::Integer, d::Integer)::Tuple
    (index ≤ 0 || index > n^d) && throw(DomainError(index, "argument must be within range [1, nᵈ]"))
    n ≤ 0 && throw(DomainError(n, "argument must be positive"))
    if d ≤ 0
        throw(DomainError(d, "argument must be postive"))
    elseif d == 1
        return (index)
    end

    remainders = Array{Int}(undef, d - 1)
    quotients = Array{Int}(undef, d)

    quotients[d] = index
    for i ∈ d-1:-1:1
        remainders[i] = ((quotients[i+1] - 1) ÷ n^i) + 1
        quotients[i] = ((quotients[i+1] - 1) % n^i) + 1
    end

    ntuple(i -> i == 1 ? quotients[i] : remainders[i-1], d)
end

# Convert a position in a d-tensor with each dimension containing n elements
# into its corresponding index. If n==3 and d==3, then a position of (2, 1, 2)
# will return 11, as (2, 1, 2) is the 11th indexed element in the tensor.
function position_index(position::Tuple, n::Integer, d::Integer)::Integer
    length(position) == d || throw(DimensionMismatch("length of the position tuple must equal d"))
    n ≤ 0 && throw(DomainError(n, "argument must be positive"))
    d ≤ 0 && throw(DomainError(d, "argument must be postive"))

    index = position[1]
    for i ∈ 2:d
        index += (position[i] - 1) * n^(i - 1)
    end

    index
end

# Given a square matrix of transition probabilities for a homogeneous Markov
# chain, produce an (k+1)-tensor listing the probabilites of each path of
# length k (i.e. the k-step transition probabilities).
function k_step_transition_probabilities(k::Integer, P::Array{Float64,2})::Array{Float64}
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
# state p₀. The path space quickly becomes large as k increases, which can
# cause memory issues with large path lengths
function exact_tvd(p₀::Integer, k::Integer, P::Array{Float64,2}, Q::Array{Float64,2})
    size(P) ≠ size(Q) && throw(DimensionMismatch("P and Q do not have the same state space"))
    size(P, 1) ≠ size(P, 2) && throw(DimensionMismatch("P and Q are not square matrices"))
    (p₀ < 0 || p₀ ≥ size(P, 1)) && throw(DomainError("p₀ should be in the range [0, √|P| - 1]"))

    n = size(P, 1)

    min_k_step_probs = min.(
        k_step_transition_probabilities(k, P),
        k_step_transition_probabilities(k, Q)
    )

    p₀_path_min_probs = min_k_step_probs[1:n:n^(k+1)]
    tvd = 1 - kahan_sum(p₀_path_min_probs)

    tvd
end
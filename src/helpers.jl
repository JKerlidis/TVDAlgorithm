# Convert the ith index of a d-tensor with each dimension containing n elements
# to a tuple representing its position in the tensor. If n==3 and d==3, then
# an index of 11 will return (2, 1, 2), because 11 is in the second row, first
# column, and at depth one in the 3-tensor.
function index_position(
    index::Integer,
    n::Integer,
    d::Integer
)::Tuple

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
function position_index(
    position::Tuple,
    n::Integer,
    d::Integer
)::Integer

    length(position) == d || throw(DimensionMismatch("length of the position tuple must equal d"))
    n ≤ 0 && throw(DomainError(n, "argument must be positive"))
    d ≤ 0 && throw(DomainError(d, "argument must be postive"))

    index = position[1]
    for i ∈ 2:d
        index += (position[i] - 1) * n^(i - 1)
    end

    index
end

# A representation of a representation that can be censored, i.e. the value may
# lie beyond the calculated range
struct CensoredObservation
    v::Float64 # Observed value
    c::Int  # Censoring flag
end

function Base.:+(
    x::CensoredObservation,
    y::CensoredObservation,
)::CensoredObservation
    CensoredObservation(x.v + y.v, x.c + y.c)
end

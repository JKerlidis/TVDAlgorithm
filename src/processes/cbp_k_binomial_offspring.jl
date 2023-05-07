# A controlled branching process (CBP) with Bin(q, m) offspring distribution,
# q ∈ ℕ₁, m ∈ [0,1], and Bin(z+K, p(z)) control function, K ∈ ℕ₁ the carrying
# capacity of the model, and p(z) = (m-1)K / ((m-2)z + mK) is a decreasing
# function of z that ensures that K is indeed a carrying capacity
struct CBPKBinomialOffspring{T<:Integer} <: TypedKBranchingProcess{T}
    K::T
    m::Integer
    q::Float64
    max_z::Integer

    function CBPKBinomialOffspring{T}(
        K::T,
        m::Integer,
        q::Float64,
        max_z::Integer
    ) where {T<:Integer}
        K ≥ 0 || throw(DomainError(K, "argument must be non-negative"))
        m ≥ 2 || throw(DomainError(m, "argument must have a value of at least 2"))
        zero(q) ≤ q ≤ one(q) || throw(DomainError(q, "argument must be in the range [0,1]"))
        max_z ≥ 0 || throw(DomainError(max_z, "argument must be non-negative"))

        new(K, m, q, max_z)
    end
end

CBPKBinomialOffspring(
    K::T,
    m::Integer,
    q::Float64,
    max_z::Integer
) where {T<:Integer} = CBPKBinomialOffspring{T}(K, m, q, max_z)

CBPKBinomialOffspring(
    K::T,
    m::Integer,
    q::Float64
) where {T<:Integer} = CBPKBinomialOffspring(K, m, q, max(trunc(Int, 3 * K), 30))

CBPKBinomialOffspring(
    K::T
) where {T<:Integer} = CBPKBinomialOffspring(K, 4, 0.25)

# Return an array of transition probabilities for the CBP
function transition_probabilities(
    d::CBPKBinomialOffspring
)::Matrix{Float64}

    p(z) = ((d.m - 1)d.K) / ((d.m - 2)z + d.m * d.K)

    P = Array{Float64}(undef, d.max_z + 1, d.max_z + 1)

    for z ∈ 0:d.max_z, c ∈ 0:d.max_z
        row_dist = Distributions.Binomial(z + d.K, p(z))

        if z == 0
            P[z+1, c+1] = c == 0 ? 1 : 0
        else
            P[z+1, c+1] = sum(
                i -> Distributions.pdf(row_dist, i) *
                     Distributions.pdf(Distributions.Binomial(d.m * i, d.q), c),
                0:(z+d.K)
            )
        end
    end

    return P
end

# Return the log-likelihood of observing a sequence of observations, under a
# given parameterisation of the CBP
function log_likelihood(
    d::CBPKBinomialOffspring,
    z_list::Vector{<:Integer},
)::Float64

    # Determine the probability ℙ((Zₙ|Zₙ₋₁=z) = j)
    p(z) = ((d.m - 1)d.K) / ((d.m - 2)z + d.m * d.K)
    b(z) = Distributions.Binomial(z + d.K, p(z))
    ℙZₙ(z, j) = z == 0 ? Int(j == 0) : sum(
        i -> Distributions.pdf(b(z), i) *
             Distributions.pdf(Distributions.Binomial(d.m * i, d.q), j),
        0:(z+d.K)
    )

    # Return the sum of negative binomial log-likelihoods, one for each
    # generation
    return sum(
        n -> log(ℙZₙ(z_list[n], z_list[n+1])),
        1:(length(z_list)-1)
    )
end

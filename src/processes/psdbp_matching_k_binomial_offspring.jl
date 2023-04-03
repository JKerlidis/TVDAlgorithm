# A population-size-dependent branching process (PSDBP) with negative binomial
# offspring distribution such that the PSDBP has a carrying capacity of K ∈ ℕ₁,
# and such that the mean and variance of the one-step distributions match that
# of the CBPKBinomialOffspring model
struct PSDBPMatchingKBinomialOffspring{T<:Integer} <: TypedBranchingProcess{T}
    K::T
    m::Integer
    q::Float64
    max_z::Integer

    function PSDBPMatchingKBinomialOffspring{T}(
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

PSDBPMatchingKBinomialOffspring(
    K::T,
    m::Integer,
    q::Float64,
    max_z::Integer
) where {T<:Integer} = PSDBPMatchingKBinomialOffspring{T}(K, m, q, max_z)

PSDBPMatchingKBinomialOffspring(
    K::T,
    m::Integer,
    q::Float64
) where {T<:Integer} = PSDBPMatchingKBinomialOffspring(K, m, q, max(trunc(Int, 3 * K), 30))

PSDBPMatchingKBinomialOffspring(
    K::T
) where {T<:Integer} = PSDBPMatchingKBinomialOffspring(K, 4, 0.25)

# Return an array of transition probabilities for the PSDBP
function transition_probabilities(
    d::PSDBPMatchingKBinomialOffspring,
)::Matrix{Float64}

    # Define the specific negative binomial distribution generating
    # Zₙ = ∑_{i=1}^z ξ(z)
    r(z, K, m) = ((z + K) * m * K) / ((m - 2) * z)
    s(z, K, m, q) = ((m - 2) * z + m * K) / ((m - 2) * z + m * K + (m - 1) * (m - 2) * q * z)
    Zₙ(z) = Distributions.NegativeBinomial(r(z, d.K, d.m), s(z, d.K, d.m, d.q))

    P = Array{Float64}(undef, d.max_z + 1, d.max_z + 1)

    for z ∈ 0:d.max_z, c ∈ 0:d.max_z
        if z == 0
            P[z+1, c+1] = c == 0 ? 1 : 0
        else
            P[z+1, c+1] = Distributions.pdf(Zₙ(z), c)
        end
    end

    return P
end

# Return the log-likelihood of observing a sequence of observations, under a
# given parameterisation of the PSDBP
function log_likelihood(
    d::PSDBPMatchingKBinomialOffspring,
    z_list::Vector{<:Integer},
)::Float64

    # Define the specific negative binomial distribution generating
    # (Zₙ|Zₙ₋₁=z) = ∑_{i=1}^z ξ(z)
    r(z, K, m) = ((z + K) * m * K) / ((m - 2) * z)
    s(z, K, m, q) = ((m - 2) * z + m * K) / ((m - 2) * z + m * K + (m - 1) * (m - 2) * q * z)
    Zₙ(z) = Distributions.NegativeBinomial(r(z, d.K, d.m), s(z, d.K, d.m, d.q))

    # Return the sum of negative binomial log-likelihoods, one for each
    # generation
    return sum(
        n -> z_list[n] == 0 ? log(z_list[n+1] == 0) : Distributions.loglikelihood(
            Zₙ(z_list[n]),
            z_list[n+1]
        ),
        1:(length(z_list)-1)
    )
end

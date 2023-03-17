# A controlled branching process (CBP) with Bin(q, m) offspring distribution,
# q ∈ ℕ₁, m ∈ [0,1], and Bin(K, p(z)) control function, K ∈ ℕ₁ the carrying
# capacity of the model, and p(z) = (m-1)K / ((m-2)z + mK) is a decreasing
# function of z that ensures that K is indeed a carrying capacity
struct CBPCarryingCapacityBinomial
    K::Integer
    m::Integer
    q::Float64
    max_z::Integer

    function CBPCarryingCapacityBinomial(
        K::Integer,
        m::Integer,
        q::Float64,
        max_z::Integer
    )
        K ≥ 0 || throw(DomainError(K, "argument must be non-negative"))
        m ≥ 2 || throw(DomainError(m, "argument must have a value of at least 2"))
        zero(q) ≤ q ≤ one(q) || throw(DomainError(q, "argument must be in the range [0,1]"))
        max_z ≥ 0 || throw(DomainError(max_z, "argument must be non-negative"))

        new(K, m, q, max_z)
    end
end

CBPCarryingCapacityBinomial(
    K::Integer,
    m::Integer,
    q::Float64
) = CBPCarryingCapacityBinomial(K, m, q, 3 * K)

CBPCarryingCapacityBinomial(
    K::Integer
) = CBPCarryingCapacityBinomial(K, 4, 0.25, 3 * K)

# Return an array of transition probabilities for the CBP
function transition_probabilities(
    d::CBPCarryingCapacityBinomial
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
    d::CBPCarryingCapacityBinomial,
    z_list::Vector{<:Integer},
)

    # Determine the probability ℙ((Zₙ|Zₙ₋₁=z) = j)
    p(z) = ((d.m - 1)d.K) / ((d.m - 2)z + d.m * d.K)
    b(z) = Distributions.Binomial(z + d.K, p(z))
    ℙZₙ(z, j) = sum(
        i -> Distributions.pdf(row_dist, i) *
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


# A population-size-dependent branching process (PSDBP) with negative binomial
# offspring distribution such that the PSDBP has a carrying capacity of K ∈ ℕ₁,
# and such that the mean and variance of the one-step distributions match that
# of the CBPCarryingCapacityBinomial model
struct PSDBPCarryingCapacityNegativeBinomial
    K::Integer
    m::Integer
    q::Float64
    max_z::Integer
    P::Matrix{Float64}

    function PSDBPCarryingCapacityNegativeBinomial(
        K::Integer,
        m::Integer,
        q::Float64,
        max_z::Integer
    )
        K ≥ 0 || throw(DomainError(K, "argument must be non-negative"))
        m ≥ 2 || throw(DomainError(m, "argument must have a value of at least 2"))
        zero(q) ≤ q ≤ one(q) || throw(DomainError(q, "argument must be in the range [0,1]"))
        max_z ≥ 0 || throw(DomainError(max_z, "argument must be non-negative"))

        new(K, m, q, max_z)
    end
end

PSDBPCarryingCapacityNegativeBinomial(
    K::Integer,
    m::Integer,
    q::Float64
) = PSDBPCarryingCapacityNegativeBinomial(K, m, q, 3 * K)

PSDBPCarryingCapacityNegativeBinomial(
    K::Integer
) = PSDBPCarryingCapacityNegativeBinomial(K, 4, 0.25, 3 * K)

# Return an array of transition probabilities for the PSDBP
function transition_probabilities(
    d::PSDBPCarryingCapacityNegativeBinomial,
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
    d::PSDBPCarryingCapacityNegativeBinomial,
    z_list::Vector{<:Integer},
)::Float64

    # Define the specific negative binomial distribution generating
    # (Zₙ|Zₙ₋₁=z) = ∑_{i=1}^z ξ(z)
    r(z, K, m) = ((z + K) * m * K) / ((m - 2) * z)
    s(z, K, m, q) = ((m - 2) * z + m * K) / ((m - 2) * z + m * K + (m - 1) * (m - 2) * q * z)
    Zₙ(z) = Distributions.NegativeBinomial(r(z, d.K, d.m), s(z, d.K, d.m, d.q)),

    # Return the sum of negative binomial log-likelihoods, one for each
    # generation
    return sum(
        n -> Distributions.loglikelihood(
            Zₙ(z_list[n]),
            z_list[n+1]
        ),
        1:(length(z_list)-1)
    )
end

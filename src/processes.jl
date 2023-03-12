# A controlled branching process (CBP) with Bin(q, m) offspring distribution,
# q ∈ ℕ₁, m ∈ [0,1], and Bin(K, p(z)) control function, K ∈ ℕ₁ the carrying
# capacity of the model, and p(z) = (m-1)K / ((m-2)z + mK) is a decreasing
# function of z that ensures that K is indeed a carrying capacity
struct CBPCarryingCapacityBinomial
    K::Integer
    m::Integer
    q::Float64
    max_z::Integer
    P::Array{Float64,2}

    function CBPCarryingCapacityBinomial(
        K::Integer,
        m::Integer,
        q::Float64,
        max_z::Integer
    )
        K ≥ 0 || throw(DomainError(K, "argument must be non-negative"))
        m ≥ 2 || throw(DomainError(m, "argument must have a value of at least 2"))
        zero(q) ≤ q ≤ one(q) || throw(DomainError(q, "argument must be in the range [0,1]"))
        max_z > 0 || throw(DomainError(max_z, "argument must be positive"))

        p(z) = ((m - 1)K) / ((m - 2)z + m * K)

        P = Array{Float64}(undef, max_z + 1, max_z + 1)

        for z ∈ 0:max_z
            for c ∈ 0:max_z
                row_dist = Distributions.Binomial(z + K, p(z))

                if z == 0
                    P[z+1, c+1] = c == 0 ? 1 : 0
                else
                    P[z+1, c+1] = sum(
                        i -> Distributions.pdf(row_dist, i) *
                             Distributions.pdf(Distributions.Binomial(m * i, q), c),
                        0:(z+K)
                    )
                end
            end
        end

        new(K, m, q, max_z, P)
    end
end

CBPCarryingCapacityBinomial(
    K::Integer
) = CBPCarryingCapacityBinomial(K, 4, 0.25, 3 * K)


# A population-size-dependent branching process (PSDBP) with negative binomial
# offspring distribution such that the PSDBP has a carrying capacity of K ∈ ℕ₁,
# and such that the mean and variance of the one-step distributions match that
# of the CBPCarryingCapacityBinomial model
struct PSDBPCarryingCapacityNegativeBinomial
    K::Integer
    m::Integer
    q::Float64
    max_z::Integer
    P::Array{Float64,2}

    function PSDBPCarryingCapacityNegativeBinomial(
        K::Integer,
        m::Integer,
        q::Float64,
        max_z::Integer
    )
        K ≥ 0 || throw(DomainError(K, "argument must be non-negative"))
        m ≥ 2 || throw(DomainError(m, "argument must have a value of at least 2"))
        zero(q) ≤ q ≤ one(q) || throw(DomainError(q, "argument must be in the range [0,1]"))
        max_z > 0 || throw(DomainError(max_z, "argument must be positive"))

        p(z) = ((m - 1)K) / ((m - 2)z + m * K)
        ξ(z) = Distributions.NegativeBinomial(
            (z + K) * m * p(z) / (z * (m * (1 - p(z)) - 1)),
            1 / (1 + q * (m - 1 - m * p(z)))
        )

        P = Array{Float64}(undef, max_z + 1, max_z + 1)

        for z ∈ 0:max_z, c ∈ 0:max_z
            if z == 0
                P[z+1, c+1] = c == 0 ? 1 : 0
            else
                P[z+1, c+1] = Distributions.pdf(ξ(z), c)
            end
        end

        new(K, m, q, max_z, P)
    end
end

PSDBPCarryingCapacityNegativeBinomial(
    K::Integer
) = PSDBPCarryingCapacityNegativeBinomial(K, 4, 0.25, 3 * K)

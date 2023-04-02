# A controlled branching process (CBP) with Poi(2) offspring distribution and
# Bin(z+M, p(z)) control function, M ∈ ℕ₁ a rate of constant immigration, and
# p(z) = K² / ((K + M)(z + K)) is a decreasing function of z that ensures that
# K > 0 is the carrying capacity the model
struct CBPKPoissonOffspring{T<:Real} <: BranchingProcess
    K::T
    M::Integer
    max_z::Integer

    function CBPKPoissonOffspring{T}(
        K::T,
        M::Integer,
        max_z::Integer
    ) where {T<:Real}
        K ≥ 0 || throw(DomainError(K, "argument must be non-negative"))
        M ≥ 0 || throw(DomainError(M, "argument must be non-negative"))
        max_z ≥ 0 || throw(DomainError(max_z, "argument must be non-negative"))

        new(K, M, max_z)
    end
end

CBPKPoissonOffspring(
    K::T,
    M::Integer,
    max_z::Integer
) where {T<:Real} = CBPKPoissonOffspring{T}(K, M, max_z)

CBPKPoissonOffspring(
    K::T,
    M::Integer
) where {T<:Real} = CBPKPoissonOffspring(K, M, max(trunc(Int, 3 * K), 30))

CBPKPoissonOffspring(
    K::T
) where {T<:Real} = CBPKPoissonOffspring(K, 2)

# Return an array of transition probabilities for the CBP
function transition_probabilities(
    d::CBPKPoissonOffspring
)::Matrix{Float64}

    p(z) = d.K^2 / ((d.K + d.M) * (z + d.K))

    P = Array{Float64}(undef, d.max_z + 1, d.max_z + 1)

    for z ∈ 0:d.max_z, c ∈ 0:d.max_z
        row_dist = Distributions.Binomial(z + d.M, p(z))

        if z == 0
            P[z+1, c+1] = c == 0 ? 1 : 0
        else
            P[z+1, c+1] = sum(
                i -> Distributions.pdf(row_dist, i) *
                     Distributions.pdf(Distributions.Poisson(2 * i), c),
                0:(z+d.M)
            )
        end
    end

    return P
end

# Return the log-likelihood of observing a sequence of observations, under a
# given parameterisation of the CBP
function log_likelihood(
    d::CBPKPoissonOffspring,
    z_list::Vector{<:Integer},
)::Float64

    # Determine the probability ℙ((Zₙ|Zₙ₋₁=z) = j)
    p(z) = d.K^2 / ((d.K + d.M) * (z + d.K))
    b(z) = Distributions.Binomial(z + d.M, p(z))
    ℙZₙ(z, j) = z == 0 ? Int(j == 0) : sum(
        i -> Distributions.pdf(b(z), i) *
             Distributions.pdf(Distributions.Poisson(2 * i), j),
        0:(z+d.M)
    )

    # Return the sum of log-likelihoods, one for each generation
    return sum(
        n -> log(ℙZₙ(z_list[n], z_list[n+1])),
        1:(length(z_list)-1)
    )
end
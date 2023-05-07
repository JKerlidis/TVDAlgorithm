# A population-size-dependent branching process (PSDBP) with negative binomial
# offspring distribution such that the PSDBP has a carrying capacity of K ∈ ℕ₁,
# and such that the mean and variance of the one-step distributions match that
# of the CBPKGeometricOffspring model. If `mean_only` is set to true, only the
# mean will be made to match
struct PSDBPMatchingKGeometricOffspring{T<:Real} <: TypedKBranchingProcess{T}
    K::T
    M::Integer
    max_z::Integer
    mean_only::Bool

    function PSDBPMatchingKGeometricOffspring{T}(
        K::T,
        M::Integer,
        max_z::Integer,
        mean_only::Bool
    ) where {T<:Real}
        K ≥ 0 || throw(DomainError(K, "argument must be non-negative"))
        M ≥ 0 || throw(DomainError(M, "argument must be non-negative"))
        max_z ≥ 0 || throw(DomainError(max_z, "argument must be non-negative"))

        new(K, M, max_z, mean_only)
    end
end

PSDBPMatchingKGeometricOffspring(
    K::T,
    M::Integer,
    max_z::Integer,
    mean_only::Bool=false
) where {T<:Real} = PSDBPMatchingKGeometricOffspring{T}(K, M, max_z, mean_only)

PSDBPMatchingKGeometricOffspring(
    K::T,
    M::Integer,
    mean_only::Bool=false
) where {T<:Real} = PSDBPMatchingKGeometricOffspring(K, M, max(trunc(Int, 3 * K), 30), mean_only)

PSDBPMatchingKGeometricOffspring(
    K::T,
    mean_only::Bool=false
) where {T<:Real} = PSDBPMatchingKGeometricOffspring(K, 2, mean_only)

# Return an array of transition probabilities for the PSDBP
function transition_probabilities(
    d::PSDBPMatchingKGeometricOffspring,
)::Matrix{Float64}

    # Define the specific negative binomial distribution generating
    # Zₙ = ∑_{i=1}^z ξ(z)
    r(z, K, M) = d.mean_only ?
                 ((z + M) * K^2) / ((K + M) * (z + K)) :
                 ((z + M) * K^2) / (2(K + M) * z + 2K * M + K^2)
    s(z, K, M) = d.mean_only ?
                 1 / 3 :
                 ((z + K) * (K + M)) / (5(K + M) * z + 5K * M + 3K^2)
    Zₙ(z) = Distributions.NegativeBinomial(r(z, d.K, d.M), s(z, d.K, d.M))

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
    d::PSDBPMatchingKGeometricOffspring,
    z_list::Vector{<:Integer},
)::Float64

    # Define the specific negative binomial distribution generating
    # (Zₙ|Zₙ₋₁=z) = ∑_{i=1}^z ξ(z)
    r(z, K, M) = d.mean_only ?
                 ((z + M) * K^2) / ((K + M) * (z + K)) :
                 ((z + M) * K^2) / (2(K + M) * z + 2K * M + K^2)
    s(z, K, M) = d.mean_only ?
                 1 / 3 :
                 ((z + K) * (K + M)) / (5(K + M) * z + 5K * M + 3K^2)
    Zₙ(z) = Distributions.NegativeBinomial(r(z, d.K, d.M), s(z, d.K, d.M))

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

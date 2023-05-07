# A population-size-dependent branching process (PSDBP) with negative binomial
# offspring distribution ξ(z) ∼ NB((z+M)/z, 1/λ), for M ∈ ℕ₁ and λ > 0, such
# that the mean and variance of the one-step distributions match that of the
# CBPSMPoissonOffspring model. The PSDBP will be a submartingale, and the mean
# number of offspring will be expected to increase each generation
struct PSDBPMatchingSMPoissonOffspring <: BranchingProcess
    M::Integer
    λ::Float64
    max_z::Integer

    function PSDBPMatchingSMPoissonOffspring(
        M::Integer,
        λ::Float64,
        max_z::Integer
    )
        M ≥ 0 || throw(DomainError(M, "argument must be non-negative"))
        λ > 0 || throw(DomainError(λ, "argument must be positive"))
        max_z ≥ 0 || throw(DomainError(max_z, "argument must be non-negative"))

        new(M, λ, max_z)
    end
end

PSDBPMatchingSMPoissonOffspring(
    max_z::Integer
) = PSDBPMatchingSMPoissonOffspring(2, 2.0, max_z)

# Return an array of transition probabilities for the PSDBP
function transition_probabilities(
    d::PSDBPMatchingSMPoissonOffspring,
)::Matrix{Float64}

    # Define the specific negative binomial distribution generating
    # Zₙ = ∑_{i=1}^z ξ(z)
    Zₙ(z) = Distributions.NegativeBinomial((z + d.M) / z, 1 / d.λ)

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
    d::PSDBPMatchingSMPoissonOffspring,
    z_list::Vector{<:Integer},
)::Float64

    # Define the specific negative binomial distribution generating
    # (Zₙ|Zₙ₋₁=z) = ∑_{i=1}^z ξ(z)
    Zₙ(z) = Distributions.NegativeBinomial((z + d.M) / z, 1 / d.λ)

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

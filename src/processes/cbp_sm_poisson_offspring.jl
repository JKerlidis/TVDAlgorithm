# A controlled branching process (CBP) with Poi(λ) offspring distribution and
# and Bin(z+M, 1/λ) control function, for M ∈ ℕ₁ a rate of constant immigration
# and λ > 0 the average number of offspring for each individual. The CBP here
# is a submartingale, and tends to grow unboundedly
struct CBPSMPoissonOffspring <: BranchingProcess
    M::Integer
    λ::Float64
    max_z::Integer

    function CBPSMPoissonOffspring(
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

CBPSMPoissonOffspring(
    max_z::Integer
) = CBPSMPoissonOffspring(2, 2.0, max_z)

# Return an array of transition probabilities for the CBP
function transition_probabilities(
    d::CBPSMPoissonOffspring
)::Matrix{Float64}

    P = Array{Float64}(undef, d.max_z + 1, d.max_z + 1)

    for z ∈ 0:d.max_z, c ∈ 0:d.max_z
        row_dist = Distributions.Binomial(z + d.M, 1 / d.λ)

        if z == 0
            P[z+1, c+1] = c == 0 ? 1 : 0
        else
            P[z+1, c+1] = sum(
                i -> Distributions.pdf(row_dist, i) *
                     Distributions.pdf(Distributions.Poisson(d.λ * i), c),
                0:(z+d.M)
            )
        end
    end

    return P
end

# Return the log-likelihood of observing a sequence of observations, under a
# given parameterisation of the CBP
function log_likelihood(
    d::CBPSMPoissonOffspring,
    z_list::Vector{<:Integer},
)::Float64

    # Determine the probability ℙ((Zₙ|Zₙ₋₁=z) = j)
    b(z) = Distributions.Binomial(z + d.M, 1 / d.λ)
    ℙZₙ(z, j) = z == 0 ? Int(j == 0) : sum(
        i -> Distributions.pdf(b(z), i) *
             Distributions.pdf(Distributions.Poisson(d.λ * i), j),
        0:(z+d.M)
    )

    # Return the sum of log-likelihoods, one for each generation
    return sum(
        n -> log(ℙZₙ(z_list[n], z_list[n+1])),
        1:(length(z_list)-1)
    )
end

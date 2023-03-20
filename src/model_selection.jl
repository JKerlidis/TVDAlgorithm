# Given a branching process, use a modification of the bisection method (under
# the assumption that this is a convex optimisation problem) to find the value
# of K in processes of same family that maximises the likelihood of observing a
# given set of data, as well as the negative log-likelihood of this maximising
# likelihood
function maximise_likelihood_in_K(
    d::BranchingProcess,
    path::Vector{Int},
    min_K::Integer=1,
    max_K::Integer=400
)::Tuple{Int,Float64}
    min_K < 0 && throw(DomainError(min_K, "min_K should be no smaller than zero"))
    max_K < min_K && throw(DomainError(max_K, "max_K should be no smaller than min_K"))

    nll_min_K = -log_likelihood(substitute_K(d, min_K), path)
    nll_max_K = -log_likelihood(substitute_K(d, max_K), path)
    while max_K > min_K
        lower_mid_K = min_K + (max_K - min_K) ÷ 2
        nll_lower_mid_K = -log_likelihood(substitute_K(d, lower_mid_K), path)

        upper_mid_K = lower_mid_K + 1
        nll_upper_mid_K = -log_likelihood(substitute_K(d, upper_mid_K), path)

        if min(nll_min_K, nll_lower_mid_K) < min(nll_upper_mid_K, nll_max_K)
            max_K = lower_mid_K
            nll_max_K = nll_lower_mid_K
        else
            min_K = upper_mid_K
            nll_min_K = nll_upper_mid_K
        end
    end

    return (min_K, nll_min_K)
end

# Given a true generating model 'd_true' and an alternative model 'd_alt',
# generate a path from the true model, and select, from the family of
# distributions taken by varying K in d_true and d_alt, the model with the
# highest likelihood of observing this path
function select_model_in_K(
    rng::Random.AbstractRNG,
    d_true::BranchingProcess,
    d_alt::BranchingProcess,
    num_steps::Integer,
    z₀::Integer=1,
    min_K::Integer=1,
    max_K::Integer=400
)::BranchingProcess

    path = Vector{Int}(undef, num_steps + 1)
    max_attempts = 10
    for i ∈ 1:max_attempts
        candidate_path = sample_path(
            rng,
            z₀,
            num_steps,
            transition_probabilities(d_true)
        )
        if !isnothing(candidate_path[num_steps+1])
            path = convert(Vector{Int}, candidate_path)
            break
        elseif i == max_attempts
            throw(ErrorException("too many attempts to generate an uncensored path"))
        end
    end

    d_true_mle_K, d_true_family_nll = maximise_likelihood_in_K(d_true, path, min_K, max_K)
    d_alt_mle_K, d_alt_family_nll = maximise_likelihood_in_K(d_alt, path, min_K, max_K)

    if d_true_family_nll < d_alt_family_nll
        return substitute_K(d_true, d_true_mle_K)
    else
        return substitute_K(d_alt, d_alt_mle_K)
    end
end

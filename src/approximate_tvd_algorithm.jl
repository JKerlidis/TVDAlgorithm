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

# Accumulator for the results of a TVD simulation
mutable struct TVDSimulationOutput
    num_trials::Integer
    num_censored::Integer
    S₁::Float64  # Cumulative sum of values
    S₂::Float64 # Cumulative sum of squared values
end

struct TVDSimulationSummary
    mean::Float64
    standard_deviation::Float64
    proportion_censored::Float64
end

function summarise(
    t::TVDSimulationOutput
)::TVDSimulationSummary

    if t.num_trials - t.num_censored == 0
        return TVDSimulationSummary(0, 0, 1)
    end

    n = t.num_trials - t.num_censored
    TVDSimulationSummary(
        t.S₁ / n,
        √(t.S₂ / n - (t.S₁ / n)^2),
        t.num_censored / t.num_trials
    )
end

# Given two homogeneous Markov processes with (possibly truncated) transition
# probabilities P and Q respectively, which are defined on the same state
# space and start from state z₀, simulate a path of length k from the process
# generated by P, and find the corresponding maximal coupling probability
# between this path and one generated from Q for each path of length 1 to k.
function simulate_coupling_probability(
    rng::Random.AbstractRNG,
    z₀::Integer,
    k::Integer,
    P::Matrix{Float64},
    Q::Matrix{Float64}
)::Array{CensoredObservation}

    size(P) ≠ size(Q) && throw(DimensionMismatch("P and Q do not have the same state space"))
    size(P, 1) ≠ size(P, 2) && throw(DimensionMismatch("P and Q are not square matrices"))
    (z₀ < 0 || z₀ ≥ size(P, 1)) && throw(DomainError(z₀, "argument should be in the range [0, √|P| - 1]"))
    k ≥ 1 || throw(DomainError(k, "path length must be at least one"))

    # Simulate k transition quantiles from P
    transition_quantiles = Random.rand(rng, k)

    # Determine the simulated path associated with these quantiles
    transition_indices = Array{Union{Int,Nothing}}(nothing, k + 1)
    transition_indices[1] = z₀
    for n ∈ 1:k
        cumulative_prob = 0.0
        for i ∈ 1:size(P, 1)
            cumulative_prob += P[transition_indices[n], i]
            if transition_quantiles[n] ≤ cumulative_prob
                transition_indices[n+1] = i
                break
            elseif i == size(P, 1)
                # If P has been truncated to before the desired quantile,
                # indicating that the data is censored, break out of both loops
                # and leave subsequent transition indices empty
                @goto path_censored
            end
        end
    end
    @label path_censored

    # Find the maximal coupling probabilities along this simulated path
    max_coupling_probs::Array{CensoredObservation} = [CensoredObservation(0, 1) for i ∈ 1:k]

    P_path_prob = 1
    Q_path_prob = 1

    for i ∈ 1:k
        if isnothing(transition_indices[i+1])
            # Return early if the path has been censored at the current
            # iteration
            return max_coupling_probs
        end

        P_path_prob *= P[transition_indices[i], transition_indices[i+1]]
        Q_path_prob *= Q[transition_indices[i], transition_indices[i+1]]

        if P_path_prob == 0
            max_coupling_probs[i] = CensoredObservation(0, 1)
        else
            max_coupling_probs[i] = CensoredObservation(max((P_path_prob - Q_path_prob) / P_path_prob, 0), 0)
        end
    end

    return max_coupling_probs
end

# Given two homogeneous Markov processes with (possibly truncated) transition
# probabilities P and Q respectively, which are defined on the same state
# space and start from state z₀, simulate the total variation distance between
# them for paths of length k, for a given number of trials. The mean
# and standard deviation are reported, as are the number and proportion of
# censored observations if the `extended_output` option is selected
function approximate_tvd(
    rng::Random.AbstractRNG,
    num_trials::Integer,
    z₀::Integer,
    k::Integer,
    P::Matrix{Float64},
    Q::Matrix{Float64}
)::Float64

    num_trials ≥ 1 || throw(DomainError(num_trials, "there must be at least one trial"))

    accumulator = CensoredObservation(0, 0)
    for i ∈ 1:num_trials
        accumulator += simulate_coupling_probability(rng, z₀, k, P, Q)[k]
    end

    accumulator.v / (num_trials - accumulator.c)
end

# Given two homogeneous Markov processes with (possibly truncated) transition
# probabilities P and Q respectively, which are defined on the same state
# space and start from state z₀, simulate the total variation distance between
# them for paths of length k, for a given number of trials. The mean
# and standard deviation are reported, as are the number and proportion of
# censored observations
function approximate_tvd_extended_output(
    rng::Random.AbstractRNG,
    num_trials::Integer,
    z₀::Integer,
    k::Integer,
    P::Array{Float64,2},
    Q::Array{Float64,2}
)::Array{TVDSimulationSummary}

    num_trials ≥ 1 || throw(DomainError(num_trials, "there must be at least one trial"))

    accumulator = [TVDSimulationOutput(num_trials, 0, 0, 0) for i ∈ 1:k]
    for i ∈ 1:num_trials
        simulated_probs = simulate_coupling_probability(rng, z₀, k, P, Q)
        for j ∈ 1:k
            accumulator[j].num_censored += simulated_probs[j].c
            accumulator[j].S₁ += simulated_probs[j].v
            accumulator[j].S₂ += simulated_probs[j].v^2
        end
    end

    summarise.(accumulator)
end

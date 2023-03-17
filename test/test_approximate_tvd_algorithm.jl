module TestApproximateTVDAlgorithm

using TVDAlgorithm
using Test
import Random
import Distributions

rng = Random.Xoshiro(1)
unif_probs = [0.25 for r ∈ 0:3, c ∈ 0:3]
bin_probs = [0.125 0.375 0.375 0.125; 0.25 0.25 0.25 0.25; 0.25 0.25 0.25 0.25; 0.125 0.375 0.375 0.125]
extinction_unif_probs = [[1 0 0 0]; [0.25 for r ∈ 0:2, c ∈ 0:3]]
extinction_bin_probs = [1 0 0 0; 0.25 0.25 0.25 0.25; 0.25 0.25 0.25 0.25; 0.125 0.375 0.375 0.125]

@testset "Test summarise function" begin
    # If every trial is censored, mean and standard error should be zero, and
    # proportion censored one
    @test summarise(TVDSimulationOutput(2, 2, -5, -5)) == TVDSimulationSummary(0, 0, 1)

    # Every trial returns the result 10 - this is deterministic so the standard
    # error should be zero
    @test summarise(TVDSimulationOutput(100, 0, 100 * 10, 100 * 10^2)) == TVDSimulationSummary(10, 0, 0)

    # Half of the trials return 0 and half 10 - the variance should be
    # 0.5 * 10² - (0.5 * 10)² = 25, and thus the standard error √(25/100) = 0.5
    @test summarise(TVDSimulationOutput(100, 0, 50 * 10, 50 * 10^2)) == TVDSimulationSummary(5, 0.5, 0)
end

@testset "Test sample_path function" begin
    struct DeterministicRNG <: Random.AbstractRNG
        sequence::Vector{Float64}
    end

    function Random.rand(dRNG::DeterministicRNG, k::Int)::Vector{Float64}
        return dRNG.sequence[1:k]
    end

    @test rand(DeterministicRNG([1, 2, 3, 4]), 3) == [1, 2, 3]

    # Test paths against deterministic quantiles
    @test sample_path(DeterministicRNG([0, 0, 0, 0]), 1, 4, unif_probs) == [1, 0, 0, 0, 0]
    @test sample_path(DeterministicRNG([0, 0.3, 0.3, 0.6]), 2, 4, unif_probs) == [2, 0, 1, 1, 2]
    @test sample_path(DeterministicRNG([0, 0.3, 0.3, 0.6]), 2, 4, bin_probs) == [2, 0, 1, 1, 2]

    # Bounds checks
    @test_nowarn sample_path(rng, 0, 4, bin_probs)
    @test_throws DomainError sample_path(rng, 5, 4, bin_probs)
end

@testset "Test simulate_coupling_probability function" begin
    @test_throws DimensionMismatch simulate_coupling_probability(
        rng, 1, 10, unif_probs, [0.25 for r ∈ 0:4, c ∈ 0:3]
    )
    for t ∈ 1:100
        simulated_prob = simulate_coupling_probability(rng, 1, 10, unif_probs, bin_probs)[10]
        @test simulated_prob.v ≥ 0 && simulated_prob.v ≤ 1
        @test simulated_prob.c == 0 || simulated_prob.c == 1
    end
end

@testset "Test approximate_tvd function" begin
    # Ensure that the TVD is symmetric
    @test round(approximate_tvd(rng, 100000, 1, 3, unif_probs, bin_probs), digits=2) ==
          round(approximate_tvd(rng, 100000, 1, 3, bin_probs, unif_probs), digits=2)

    # Sense checks on the TVD bound
    @test approximate_tvd(rng, 100, 0, 5, extinction_unif_probs, extinction_bin_probs) == 0
    @test approximate_tvd(rng, 100, 1, 5, extinction_unif_probs, extinction_bin_probs) > 0
    @test round(approximate_tvd(rng, 10000, 1, 5, extinction_unif_probs, extinction_bin_probs), digits=2) ==
          round(approximate_tvd(rng, 10000, 2, 5, extinction_unif_probs, extinction_bin_probs), digits=2)
    @test approximate_tvd(rng, 100, 1, 5, extinction_unif_probs, extinction_bin_probs) <
          approximate_tvd(rng, 100, 3, 5, extinction_unif_probs, extinction_bin_probs)
    @test_throws DomainError approximate_tvd(rng, 100, 4, 5, extinction_unif_probs, extinction_bin_probs)

    @test round(approximate_tvd(rng, 10000000, 1, 3, unif_probs, bin_probs), digits=3) == 0.188
    @test round(approximate_tvd(rng, 10000000, 1, 3, extinction_bin_probs, extinction_unif_probs), digits=3) == 0.094
    @test round(approximate_tvd(rng, 10000000, 0, 3, extinction_bin_probs, unif_probs), digits=3) == 0.984
end

end

module TestApproximateTVDAlgorithm

using TVDAlgorithm
using Test
import Random
import Distributions

rng = Random.Xoshiro(1)
unif_probs = [0.25 for r ∈ 0:3, c ∈ 0:3]
bin_probs = [0.125 0.375 0.375 0.125; 0.25 0.25 0.25 0.25; 0.25 0.25 0.25 0.25; 0.125 0.375 0.375 0.125]

@testset "Test simulate_coupling_probability function" begin
    @test_throws DimensionMismatch simulate_coupling_probability(
        rng, 1, 10, unif_probs, [0.25 for r ∈ 0:4, c ∈ 0:3]
    )
    for t ∈ 1:50
        simulated_prob = simulate_coupling_probability(rng, 1, 10, unif_probs, bin_probs)[10]
        @test simulated_prob.v ≥ 0 && simulated_prob.v ≤ 1
        @test simulated_prob.c == 0 || simulated_prob.c == 1
    end
end

@testset "Test approximate_tvd function" begin
    @test round(approximate_tvd(rng, 10000000, 1, 3, unif_probs, bin_probs), digits=3) ==
          round(approximate_tvd(rng, 10000000, 1, 3, bin_probs, unif_probs), digits=3)
    @test round(approximate_tvd(rng, 10000000, 1, 3, unif_probs, bin_probs), digits=3) == 0.281
end

end

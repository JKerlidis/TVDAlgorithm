module TestProcesses

using TVDAlgorithm
using Test
import Distributions

cbp_3 = CBPCarryingCapacityBinomialOffspring(3, 10, 0.1, 5)
psdbp_3 = PSDBPMatchingKBinomialOffspring(3, 10, 0.1, 5)
cbp_10 = CBPCarryingCapacityBinomialOffspring(10, 5, 0.2)
psdbp_10 = PSDBPMatchingKBinomialOffspring(10, 5, 0.2)

@testset "Test substitute_K function" begin
    @test substitute_K(cbp_3, 10) == CBPCarryingCapacityBinomialOffspring(10, 10, 0.1, 5)
    @test substitute_K(psdbp_3, 10) == PSDBPMatchingKBinomialOffspring(10, 10, 0.1, 5)
end

@testset "Test transition_probabilities function" begin
    # Sense checks
    @test transition_probabilities(cbp_3)[1, :] == [1, 0, 0, 0, 0, 0]
    @test transition_probabilities(psdbp_3)[1, :] == [1, 0, 0, 0, 0, 0]
    @test all(sum(r) ≤ 1 for r ∈ eachrow(transition_probabilities(cbp_10)))
    @test all(sum(r) ≤ 1 for r ∈ eachrow(transition_probabilities(psdbp_10)))

    # Specific calculated values
    @test round(transition_probabilities(cbp_3)[2, 3], digits=3) == 0.217
    @test round(transition_probabilities(psdbp_3)[4, 1], digits=3) == 0.080
end

@testset "Test log_likelihood function" begin
    # Sense checks
    @test log_likelihood(cbp_3, [4, 2, 1, 3, 0]) < 0
    @test log_likelihood(psdbp_3, [4, 2, 1, 3, 0]) < 0
    @test log_likelihood(cbp_3, [0, 0]) == 0
    @test log_likelihood(cbp_3, [0, 1]) == -Inf
    @test log_likelihood(psdbp_3, [0, 0]) == 0
    @test log_likelihood(psdbp_3, [0, 1]) == -Inf

    # Specific calculated values (to match the transition_probabilities checks)
    @test round(log_likelihood(cbp_3, [1, 2]), digits=3) == round(log(0.2165722), digits=3)
    @test round(log_likelihood(psdbp_3, [3, 0]), digits=3) == round(log(0.0801751), digits=3)
end

end

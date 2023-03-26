module TestProcesses

using TVDAlgorithm
using Test
import Distributions

cbp_bin = CBPKBinomialOffspring(3, 10, 0.1, 5)
psdbp_bin = PSDBPMatchingKBinomialOffspring(3, 10, 0.1, 5)
cbp_geom = CBPKGeometricOffspring(10, 2, 5)
psdbp_geom = PSDBPMatchingKGeometricOffspring(10, 2, 5)

@testset "Test substitute_K function" begin
    @test substitute_K(cbp_bin, 10) == CBPKBinomialOffspring(10, 10, 0.1, 5)
    @test substitute_K(psdbp_bin, 10) == PSDBPMatchingKBinomialOffspring(10, 10, 0.1, 5)
    @test substitute_K(cbp_geom, 6) == CBPKGeometricOffspring(6, 2, 5)
    @test substitute_K(psdbp_geom, 6, true) == PSDBPMatchingKGeometricOffspring(6, 2, 30)
end

@testset "Test transition_probabilities function" begin
    # Sense checks
    @test transition_probabilities(cbp_bin)[1, :] == [1, 0, 0, 0, 0, 0]
    @test transition_probabilities(psdbp_bin)[1, :] == [1, 0, 0, 0, 0, 0]
    @test transition_probabilities(cbp_geom)[1, :] == [1, 0, 0, 0, 0, 0]
    @test transition_probabilities(psdbp_geom)[1, :] == [1, 0, 0, 0, 0, 0]
    @test all(sum(r) ≤ 1 for r ∈ eachrow(transition_probabilities(cbp_bin)))
    @test all(sum(r) ≤ 1 for r ∈ eachrow(transition_probabilities(psdbp_bin)))
    @test all(sum(r) ≤ 1 for r ∈ eachrow(transition_probabilities(cbp_geom)))
    @test all(sum(r) ≤ 1 for r ∈ eachrow(transition_probabilities(psdbp_geom)))

    # Specific calculated values
    @test round(transition_probabilities(cbp_bin)[2, 3], digits=3) == 0.217
    @test round(transition_probabilities(psdbp_bin)[4, 1], digits=3) == 0.080
    @test round(transition_probabilities(cbp_geom)[2, 1], digits=3) == 0.121
    @test round(transition_probabilities(psdbp_geom)[4, 4], digits=3) == 0.102
end

@testset "Test log_likelihood function" begin
    # Sense checks
    @test log_likelihood(cbp_bin, [4, 2, 1, 3, 0]) < 0
    @test log_likelihood(psdbp_bin, [4, 2, 1, 3, 0]) < 0
    @test log_likelihood(cbp_geom, [4, 2, 1, 3, 0]) < 0
    @test log_likelihood(psdbp_geom, [4, 2, 1, 3, 0]) < 0
    @test log_likelihood(cbp_bin, [0, 0]) == 0
    @test log_likelihood(cbp_bin, [0, 1]) == -Inf
    @test log_likelihood(psdbp_bin, [0, 0]) == 0
    @test log_likelihood(psdbp_bin, [0, 1]) == -Inf
    @test log_likelihood(cbp_geom, [0, 0]) == 0
    @test log_likelihood(cbp_geom, [0, 1]) == -Inf
    @test log_likelihood(psdbp_geom, [0, 0]) == 0
    @test log_likelihood(psdbp_geom, [0, 1]) == -Inf

    # Specific calculated values (to match the transition_probabilities checks)
    @test round(log_likelihood(cbp_bin, [1, 2]), digits=5) == round(log(0.2165722), digits=5)
    @test round(log_likelihood(psdbp_bin, [3, 0]), digits=5) == round(log(0.0801751), digits=5)
    @test round(log_likelihood(cbp_geom, [1, 0]), digits=5) == round(log(0.1212503), digits=5)
    @test round(log_likelihood(psdbp_geom, [3, 3]), digits=5) == round(log(0.1015597), digits=5)
end

end

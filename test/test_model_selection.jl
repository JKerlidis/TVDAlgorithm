module TestModelSelection

using TVDAlgorithm
using Random
using Test

cbp_bin = CBPKBinomialOffspring(10, 5, 0.2)
psdbp_bin = PSDBPMatchingKBinomialOffspring(10, 5, 0.2)
path_bin = [1, 8, 12, 6, 11, 5, 10, 8, 7, 6, 7, 7, 9, 6, 11, 14, 8, 5, 15, 8, 12, 9, 6, 8, 6, 9, 13, 17, 18, 15, 10]

cbp_poi = CBPKPoissonOffspring(10.0, 2)
psdbp_poi = PSDBPMatchingKPoissonOffspring(10.0, 2)
path_poi = [1, 5, 6, 8, 8, 18, 6, 7, 6, 17, 11, 12, 16, 5, 14, 15, 16, 14, 7, 9, 9, 3, 5, 7, 15, 12, 16, 10, 11, 15, 14]

@testset "Test maximise_likelihood_in_K function" begin
    # Sense checks
    @test_throws DomainError maximise_likelihood_in_K(cbp_bin, path_bin, 3, 2)
    @test_throws DomainError maximise_likelihood_in_K(cbp_poi, path_poi, 3, 2)
    @test maximise_likelihood_in_K(psdbp_bin, path_bin, 5, 5) ==
          (5, -log_likelihood(substitute_K(psdbp_bin, 5), path_bin))
    @test maximise_likelihood_in_K(psdbp_poi, path_poi, 5, 5) ==
          (5.0, -log_likelihood(substitute_K(psdbp_poi, 5), path_poi))

    # Specific values
    @test maximise_likelihood_in_K(cbp_bin, path_bin)[1] == 10
    @test maximise_likelihood_in_K(psdbp_bin, path_bin)[1] == 10

    cbp_poi_mle, cbp_poi_nll = maximise_likelihood_in_K(cbp_poi, path_poi)
    @test cbp_poi_nll < -log_likelihood(substitute_K(cbp_poi, floor(Int, cbp_poi_mle)), path_poi)
    @test cbp_poi_nll < -log_likelihood(substitute_K(cbp_poi, ceil(Int, cbp_poi_mle)), path_poi)

    psdbp_poi_mle, psdbp_poi_nll = maximise_likelihood_in_K(psdbp_poi, path_poi)
    @test psdbp_poi_nll < -log_likelihood(substitute_K(psdbp_poi, floor(Int, psdbp_poi_mle)), path_poi)
    @test psdbp_poi_nll < -log_likelihood(substitute_K(psdbp_poi, ceil(Int, psdbp_poi_mle)), path_poi)
end

@testset "Test select_model_in_K function" begin
    # Make sure integer-valued branching processes are routed to the integer
    # version of `maximise_likelihood_in_K`, and analagously for floating-point
    # valued branching processes
    @test select_model_in_K(Xoshiro(1), transition_probabilities(psdbp_bin), psdbp_bin, cbp_bin, 30) isa
          TypedBranchingProcess{Int}
    @test select_model_in_K(Xoshiro(1), transition_probabilities(psdbp_poi), psdbp_poi, cbp_poi, 30) isa
          TypedBranchingProcess{Float64}
end

end

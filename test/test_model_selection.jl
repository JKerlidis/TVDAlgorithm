module TestModelSelection

using TVDAlgorithm
using Test

cbp_bin = CBPKBinomialOffspring(10, 5, 0.2)
psdbp_bin = PSDBPMatchingKBinomialOffspring(10, 5, 0.2)
path = [1, 8, 12, 6, 11, 5, 10, 8, 7, 6, 7, 7, 9, 6, 11, 14, 8, 5, 15, 8, 12, 9, 6, 8, 6, 9, 13, 17, 18, 15, 10]

@testset "Test maximise_likelihood_in_K function" begin
    # Sense checks
    @test_throws DomainError maximise_likelihood_in_K(cbp_bin, path, 3, 2)
    @test maximise_likelihood_in_K(psdbp_bin, path, 5, 5) == (5, -log_likelihood(substitute_K(psdbp_bin, 5), path))

    # Specific values
    @test maximise_likelihood_in_K(cbp_bin, path)[1] == 10
    @test maximise_likelihood_in_K(psdbp_bin, path)[1] == 10
end

end

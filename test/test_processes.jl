module TestProcesses

using TVDAlgorithm
using Test
import Distributions

cbp_3 = CBPCarryingCapacityBinomial(3, 10, 0.1, 5)
psdbp_3 = PSDBPCarryingCapacityNegativeBinomial(3, 10, 0.1, 5)
cbp_10 = CBPCarryingCapacityBinomial(10, 5, 0.2)
psdbp_10 = PSDBPCarryingCapacityNegativeBinomial(10, 5, 0.2)

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

end

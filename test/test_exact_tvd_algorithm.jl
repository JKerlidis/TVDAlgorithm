module TestExactTVDAlgorithm

using TVDAlgorithm
using Test
import Distributions

@testset "Test exact_tvd function" begin
    unif_probs = [0.25 for r ∈ 0:3, c ∈ 0:3]
    bin_probs = [0.125 0.375 0.375 0.125; 0.25 0.25 0.25 0.25; 0.25 0.25 0.25 0.25; 0.125 0.375 0.375 0.125]

    @test exact_tvd(1, 3, bin_probs, unif_probs) ≈ exact_tvd(1, 3, unif_probs, bin_probs)
    @test round(exact_tvd(1, 3, bin_probs, unif_probs), digits=3) == 0.281
end

end

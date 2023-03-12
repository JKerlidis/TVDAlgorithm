module TestExactTVDAlgorithm

using TVDAlgorithm
using Test
import Distributions

@testset "Test exact_tvd function" begin
    unif_probs = [0.25 for r ∈ 0:3, c ∈ 0:3]
    bin_probs = [Distributions.pdf(Distributions.Binomial(3, 0.5), c) for r ∈ 0:3, c ∈ 0:3]

    @test exact_tvd(1, 3, bin_probs, unif_probs) ≈ exact_tvd(1, 3, unif_probs, bin_probs)
    @test exact_tvd(1, 3, bin_probs, unif_probs) ≈ 0.34375
end

end

module TestExactTVDAlgorithm

using TVDAlgorithm
using Test
import Distributions

unif_probs = [0.25 for r ∈ 0:3, c ∈ 0:3]
bin_probs = [0.125 0.375 0.375 0.125; 0.25 0.25 0.25 0.25; 0.25 0.25 0.25 0.25; 0.125 0.375 0.375 0.125]
extinction_unif_probs = [[1 0 0 0]; [0.25 for r ∈ 0:2, c ∈ 0:3]]
extinction_bin_probs = [1 0 0 0; 0.25 0.25 0.25 0.25; 0.25 0.25 0.25 0.25; 0.125 0.375 0.375 0.125]

@testset "Test exact_tvd function" begin
      # Ensure that the TVD is symmetric
      @test exact_tvd(1, 3, bin_probs, unif_probs) ≈ exact_tvd(1, 3, unif_probs, bin_probs)

      # Sense checks on the TVD bound
      @test exact_tvd(0, 5, extinction_unif_probs, extinction_bin_probs) == 0
      @test exact_tvd(1, 5, extinction_unif_probs, extinction_bin_probs) > 0
      @test exact_tvd(1, 5, extinction_unif_probs, extinction_bin_probs) ==
            exact_tvd(2, 5, extinction_unif_probs, extinction_bin_probs)
      @test exact_tvd(1, 5, extinction_unif_probs, extinction_bin_probs) <
            exact_tvd(3, 5, extinction_unif_probs, extinction_bin_probs)
      @test_throws DomainError exact_tvd(4, 5, extinction_unif_probs, extinction_bin_probs)

      # Ensure that exact TVD values match the approximate TVD calculation
      @test round(exact_tvd(1, 3, bin_probs, unif_probs), digits=3) == 0.188
      @test round(exact_tvd(1, 3, extinction_bin_probs, extinction_unif_probs), digits=3) == 0.094
      @test round(exact_tvd(0, 3, extinction_bin_probs, unif_probs), digits=3) == 0.984
end

end

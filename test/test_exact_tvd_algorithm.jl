module TestExactTVDAlgorithm

using TVDAlgorithm
using Test
import Distributions

@testset "Test index_position function" begin
    @test_throws DomainError index_position(20, 2, 4)
    @test index_position(11, 3, 3) == (2, 1, 2)
    @test index_position(18, 3, 3) == (3, 3, 2)
    @test index_position(56, 3, 4) == (2, 1, 1, 3)
end

@testset "Test position_index function" begin
    @test_throws DimensionMismatch position_index((1, 2, 3), 3, 4)
    @test position_index((2, 1, 2), 3, 3) == 11
    @test position_index((3, 3, 2), 3, 3) == 18
    @test position_index((2, 1, 1, 3), 3, 4) == 56
end

@testset "Test kahan_sum function" begin
    @test kahan_sum([1, 2, 3, 4]) == 10
    @test kahan_sum(repeat([1 / 700, 3 / 700], 2100)) == 12.0
end

@testset "Test exact_tvd function" begin
    unif_probs = [0.25 for r ∈ 0:3, c ∈ 0:3]
    bin_probs = [Distributions.pdf(Distributions.Binomial(3, 0.5), c) for r ∈ 0:3, c ∈ 0:3]

    @test exact_tvd(1, 3, bin_probs, unif_probs) ≈ exact_tvd(1, 3, unif_probs, bin_probs)
    @test exact_tvd(1, 3, bin_probs, unif_probs) ≈ 0.34375
end

end

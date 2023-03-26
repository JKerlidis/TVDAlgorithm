module TestExactTVDAlgorithm

using TVDAlgorithm
using Test

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

@testset "Test CensoredObservation" begin
    @test CensoredObservation(3, 1) + CensoredObservation(5, 0) == CensoredObservation(8, 1)
end

end

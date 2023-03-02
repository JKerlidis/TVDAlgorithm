module TestCompensatedSummation

using TVDAlgorithm
using Test

@testset "Test kahan_sum function" begin
    @test kahan_sum([1, 2, 3, 4]) == 10
    @test kahan_sum(repeat([1 / 700, 3 / 700], 2100)) == 12.0
end

end

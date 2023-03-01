using TVDAlgorithm
using Test

const tests = [
    "test_exact_tvd_algorithm"
]

@testset "TVDAlgorithm" begin
    @testset "Test $t" for t ∈ tests
        include("$t.jl")
    end
end

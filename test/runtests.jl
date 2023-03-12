using TVDAlgorithm
using Test

const tests = [
    "test_helpers",
    "test_exact_tvd_algorithm",
    "test_approximate_tvd_algorithm"
]

@testset "TVDAlgorithm" begin
    @testset "Test $t" for t ∈ tests
        include("$t.jl")
    end
end

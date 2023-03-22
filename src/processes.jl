abstract type BranchingProcess end

const processes = [
    "cbp_k_binomial_offspring",
    "cbp_k_geometric_offspring",
    "psdbp_matching_k_binomial_offspring",
    "psdbp_matching_k_geometric_offspring"
]

function substitute_K(
    d::T,
    new_K::Integer
)::T where {T<:BranchingProcess}
    return T(
        [
            key == :K ? new_K : getfield(d, key) for key ∈ fieldnames(T)
        ]...
    )
end

for p ∈ processes
    include(joinpath("processes", "$p.jl"))
end

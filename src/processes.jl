abstract type BranchingProcess end

const processes = [
    "cbp_k_binomial_offspring",
    "cbp_k_geometric_offspring",
    "psdbp_matching_k_binomial_offspring",
    "psdbp_matching_k_geometric_offspring"
]

# Return a copy of the branching process with a new carrying capacity (K) set.
# Optionally also reset the process' value of `max_z` to reflect this new
# carrying capacity
function substitute_K(
    d::T,
    new_K::Integer,
    reset_max_z::Bool=false
)::T where {T<:BranchingProcess}

    return T(
        [
            if key == :K
                new_K
            elseif key == :max_z
                reset_max_z ? max(3 * new_K, 30) : getfield(d, key)
            else
                getfield(d, key)
            end
            for key ∈ fieldnames(T)
        ]...
    )
end

for p ∈ processes
    include(joinpath("processes", "$p.jl"))
end

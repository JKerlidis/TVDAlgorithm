abstract type BranchingProcess end

const processes = [
    "cbp_k_binomial_offspring",
    "cbp_k_geometric_offspring",
    "cbp_k_poisson_offspring",
    "psdbp_matching_k_binomial_offspring",
    "psdbp_matching_k_geometric_offspring",
    "psdbp_matching_k_poisson_offspring"
]

# Return a copy of the branching process with a new carrying capacity (K) set.
# Optionally also reset the process' value of `max_z` to reflect this new
# carrying capacity
function substitute_K(
    d::BranchingProcess,
    new_K::Real,
    reset_max_z::Bool=false
)::BranchingProcess

    T = typeof(d)

    return T(
        [
            if key == :K
                getfield(d, key) isa Integer ? trunc(Int, new_K) : convert(typeof(getfield(d, key)), new_K)
            elseif key == :max_z
                reset_max_z ? max(trunc(Int, 3 * new_K), 30) : getfield(d, key)
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

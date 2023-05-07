abstract type BranchingProcess end
abstract type TypedKBranchingProcess{T} <: BranchingProcess end

const processes = [
    "cbp_k_binomial_offspring",
    "cbp_k_geometric_offspring",
    "cbp_k_poisson_offspring",
    "cbp_sm_poisson_offspring",
    "psdbp_matching_k_binomial_offspring",
    "psdbp_matching_k_geometric_offspring",
    "psdbp_matching_k_poisson_offspring",
    "psdbp_matching_sm_poisson_offspring",
]

# Return a copy of the branching process with a new carrying capacity (K) set.
# Optionally also reset the process' value of `max_z` to reflect this new
# carrying capacity
function substitute_K(
    d::TypedKBranchingProcess{T},
    new_K::Real,
    reset_max_z::Bool=false
) where {T<:Real}

    if T <: Integer
        set_K = trunc(T, new_K)
        set_max_z = max(3 * set_K, 30)
    else
        set_K = convert(T, new_K)
        set_max_z = max(trunc(Int, 3 * new_K), 30)
    end

    S = typeof(d)
    fields = [
        if key == :K
            set_K
        elseif key == :max_z
            reset_max_z ? set_max_z : getfield(d, key)
        else
            getfield(d, key)
        end
        for key ∈ fieldnames(S)
    ]

    return S(fields...)
end

for p ∈ processes
    include(joinpath("processes", "$p.jl"))
end

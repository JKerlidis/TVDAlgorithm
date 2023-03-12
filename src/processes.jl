struct CBPBevertonHoltBinomial
    K::Integer
    m::Integer
    q::Float64
    max_z::Integer
    P::Array{Float64,2}
    function CBPBevertonHoltBinomial(K::Integer, m::Integer, q::Float64, max_z::Integer)
        K ≥ 0 || throw(DomainError(K, "argument must be non-negative"))
        m = max(m, 0)
        zero(q) ≤ q ≤ one(q) || throw(DomainError(q, "argument must be in the range [0,1]"))
        max_z ≥ 0 || throw(DomainError(max_z, "argument must be non-negative"))

        p(z) = (3K) / (2z + 4K)

        P = zeros(Float64, max_z + 1, max_z + 1)


        new(K, m, q, max_z)
    end
end
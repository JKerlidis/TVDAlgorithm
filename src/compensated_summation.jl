# An implementation of the Kahan compensated summation algorithm, to minimise
# the impact of floating-point-error when summing many floats
function kahan_sum(summands::Array{T})::T where {T}
    sum = 0.0
    compensator = 0.0

    for s ∈ summands
        t = sum + s
        if abs(sum) ≥ abs(s)
            compensator += (sum - t) + s
        else
            compensator += (s - t) + sum
        end
        sum = t
    end

    sum + compensator
end

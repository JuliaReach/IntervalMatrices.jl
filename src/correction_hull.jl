"""
    correction_hull(A::IntervalMatrix{T}, t, p) where {T}

Compute the correction term for the convex hull of a point and its linear map
with an interval matrix in order to contain all trajectories of a linear system.

### Input

- `A` -- interval matrix
- `t` -- non-negative time value
- `p` -- order of the approximation

### Output

An interval matrix representing the correction term.

### Algorithm

See Theorem 3 in [1].

[1] M. Althoff, O. Stursberg, M. Buss. Reachability Analysis of Linear Systems
with Uncertain Parameters and Inputs. CDC 2007.
"""
function correction_hull(A::IntervalMatrix{T}, t, p) where {T}
    # initialize interval matrix with zero intervals
    m, n = size(A)
    F = IntervalMatrix(zeros(Interval{T}, m, n))

    A2i = A
    fac_i = 1
    t2i = t
    @inbounds for i in 2:p
        t2i *= t
        left = (one(T) / i^(i/i-1) - one(T) / i^(1/i-1)) * t2i
        itv = Interval(left, zero(T))
        A2i = A2i * A
        fac_i *= i
        F += itv * A2i * (1/fac_i)
    end
    F += _expm_remainder(A, t, p)
    return F
end

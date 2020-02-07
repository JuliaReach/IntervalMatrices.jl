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
    F = _exp_remainder(A, t, p)
    _correction_loop!(F, A, A, t, p)
    return F
end

"""
    input_correction(A::IntervalMatrix{T}, t, p) where {T}

Compute the *input correction matrix* for discretizing an inhomogeneous affine
dynamical system with an interval matrix and an input domain not containing the
origin.

### Input

- `A` -- interval matrix
- `t` -- non-negative time value
- `p` -- order of the Taylor approximation

### Output

An interval matrix representing the correction matrix.

### Algorithm

See Proposition 3.4 in [1].

[1] M. Althoff. Reachability analysis and its application to the safety
assessment of automonous cars. 2010.
"""
function input_correction(A::IntervalMatrix{T}, t, p) where {T}
    n = checksquare(A)

    E = _exp_remainder(A, t, p; n=n)
    F = E / opnorm(A, Inf)

    id = IntervalMatrix(Diagonal(fill(Interval(one(T)), n)))  # identity matrix
    _correction_loop!(F, A, id, t, p)

    return F
end

# inputs:
# F:  interval matrix sum (modified)
# A:  interval matrix
# Aⁱ: interval matrix factor in iteration i to start with; examples:
#     if Aⁱ == A, then this function will multiply with Aⁱ
#     if Aⁱ == id, then this function will multiply with Aⁱ⁻¹
# t:  time step
# p:  Taylor-approximation order
# outputs:
# F is modified by adding \\sum_{i=2}^p [...., 0] * t^i * (Aⁱ * A^i)
function _correction_loop!(F, A::IntervalMatrix{T}, Aⁱ, t, p) where {T}
    i! = 1
    tⁱ = t
    @inbounds for i in 2:p
        tⁱ *= t
        left = (one(T) / i^(i/(i-1)) - one(T) / i^(1/(i-1))) * tⁱ
        itv = Interval(left, zero(T))
        Aⁱ *= A
        i! *= i
        F += itv * Aⁱ / i!
    end
    return F
end

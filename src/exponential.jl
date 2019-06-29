"""
    quadratic_expansion(A::IntervalMatrix, t, p)

Exactly compute the quadratic formula ``At + \\frac{1}{2}A^2t^2`` using interval
arithmetics.

### Input

- `A` -- interval matrix
- `t` -- non-negative time value

### Algorithm

See Lemma 1 in *Reachability Analysis of Linear Systems with Uncertain Parameters
and Inputs* by M. Althoff, O. Stursberg, M. Buss.
"""
function quadratic_expansion(A::IntervalMatrix, t)
    n = LinearAlgebra.checksquare(A)

    @inline function κ(aii, t)
        if -1/t ∈ aii
            return -1/2
        else
            return min(aii.left*t+aii.left^2*t^2/2, aii.right*t+aii.right^2*t^2/2)
        end
    end

    W = similar(A)

    for i in 1:n
        for j in 1:n
            if i ≠ j
                W[i, j] = A[i, j] * (t + (A[i, i] + A[j, j])*(t^2/2.))
                S = 0
                for k in 1:n
                    if k ≠ i && k ≠ j
                        S = S + A[i, k] * A[k, j]
                    end
                end
                W[i, j] = W[i, j] + S * (t^2/2.)
            else
                u = A[i, i].left * t + A[i, i].left^2 * t^2/2
                v = A[i, i].right * t + A[i, i].right^2 * t^2/2
                W[i, i] = ClosedInterval(κ(A[i, i], t), max(u, v))
                S = 0
                for k in 1:n
                    if k ≠ i
                        S = S + A[i, k] * A[k, i]
                    end
                end
                W[i, i] = W[i, i] + S * (t^2/2.)
            end
        end
    end
    return W
end

"""
    expm_overapproximation(M::IntervalMatrix{T, <: AbstractInterval{T}}, t, p) where {T}

Overapproximation of the exponential of an interval matrix.

### Input

- `A` -- interval matrix
- `t` -- non-negative time value
- `p` -- order of the appproximation

### Algorithm

See Theorem 1 in *Reachability Analysis of Linear Systems with Uncertain
Parameters and Inputs* by M. Althoff, O. Stursberg, M. Buss.
"""
function expm_overapproximation(A::IntervalMatrix{T, <: AbstractInterval{T}}, t, p) where {T}
    nA = opnorm(A, Inf)
    c = nA * t / (p + 2)
    @assert c < 1

    n = size(A, 1)
    Γ = IntervalMatrix(fill(zero(T)±one(T), (n , n)))
    E = Γ * ((nA*t)^(p+1) * (1/factorial(p + 1) * 1/(1-c)))

    S = IntervalMatrix(fill(zero(T)±zero(T), (n , n)))
    Ai = A * A
    fact_num = t^2
    fact_denom = 2
    for i in 3:p
        fact_num *= t
        fact_denom *= i
        Ai = Ai * A
        S = S + Ai * (fact_num / fact_denom)
    end
    W = quadratic_expansion(A, t)
    res = W + S + E

    # add identity matrix implicitly
    for i in 1:n
        res[i, i] += one(T)
    end
    return res
end

"""
    expm_underapproximation(M::IntervalMatrix{T, <: AbstractInterval{T}}, t, p) where {T}

Overapproximation of the exponential of an interval matrix.

### Input

- `A` -- interval matrix
- `t` -- non-negative time value
- `p` -- order of the appproximation

### Algorithm

See Theorem 2 in *Reachability Analysis of Linear Systems with Uncertain
Parameters and Inputs* by M. Althoff, O. Stursberg, M. Buss.
"""
function expm_underapproximation(A::IntervalMatrix{T, <: AbstractInterval{T}}, t, p) where {T}
    n = size(A, 1)

    Y = zeros(n, n)
    LA = left(A)
    Ail = LA * LA
    Z = zeros(n, n)
    RA = right(A)
    Air = RA * RA
    fact_num = t^2
    fact_denom = 2
    for i in 3:p
        fact_num *= t
        fact_denom *= i
        fact = fact_num / fact_denom
        Ail = Ail * LA
        Y = Y + Ail * fact
        Air = Air * RA
        Z = Z + Air * fact
    end

    B = IntervalMatrix(Matrix{Interval{:closed, :closed, Float64}}(undef, n , n))
    for j in 1:n
        for i in 1:n
            B[i, j] = ClosedInterval(min(Y[i, j], Z[i, j]), max(Y[i, j], Z[i, j]))
        end
    end

    W = quadratic_expansion(A, t)
    res = W + B

    # add identity matrix implicitly
    for i in 1:n
        res[i, i] += one(T)
    end
    return res
end

using LinearAlgebra: checksquare

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
    n = checksquare(A)

    t2d2 = t^2/2
    @inline function κ(aii, t)
        if -1/t ∈ aii
            return -1/2
        else
            a = inf(aii) * t + inf(aii)^2 * t2d2
            b = sup(aii) * t + sup(aii)^2 * t2d2
            return min(a, b)
        end
    end

    W = similar(A)

    @inbounds for j in 1:n
        for i in 1:n
            S = 0
            if i ≠ j
                for k in 1:n
                    if k ≠ i && k ≠ j
                        S = S + A[i, k] * A[k, j]
                    end
                end
                W[i, j] = A[i, j] * (t + (A[i, i] + A[j, j]) * t2d2) + S * t2d2
            else
                for k in 1:n
                    if k ≠ i
                        S = S + A[i, k] * A[k, j]
                    end
                end
                u = inf(A[i, i]) * t + inf(A[i, i])^2 * t2d2
                v = sup(A[i, i]) * t + sup(A[i, i])^2 * t2d2
                W[i, i] = Interval(κ(A[i, i], t), max(u, v)) + S * t2d2
            end
        end
    end
    return W
end

"""
    expm_overapproximation(M::IntervalMatrix{T, Interval{T}}, t, p) where {T}

Overapproximation of the exponential of an interval matrix.

### Input

- `A` -- interval matrix
- `t` -- non-negative time value
- `p` -- order of the approximation

### Algorithm

See Theorem 1 in *Reachability Analysis of Linear Systems with Uncertain
Parameters and Inputs* by M. Althoff, O. Stursberg, M. Buss.
"""
function expm_overapproximation(A::IntervalMatrix{T, Interval{T}}, t, p) where {T}
    n = checksquare(A)

    E = _expm_remainder(A, t, p; n=n)
    S = IntervalMatrix(zeros(Interval{T}, n, n))
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

function _expm_remainder(A::IntervalMatrix{T}, t, p; n=checksquare(A)) where {T}
    nA = opnorm(A, Inf)
    c = nA * t / (p + 2)
    @assert c < 1
    Γ = IntervalMatrix(fill(zero(T)±one(T), (n , n)))
    return Γ * ((nA*t)^(p+1) * (1/factorial(p + 1) * 1/(1-c)))
end

"""
    expm_underapproximation(M::IntervalMatrix{T, Interval{T}}, t, p) where {T}

Overapproximation of the exponential of an interval matrix.

### Input

- `A` -- interval matrix
- `t` -- non-negative time value
- `p` -- order of the approximation

### Algorithm

See Theorem 2 in *Reachability Analysis of Linear Systems with Uncertain
Parameters and Inputs* by M. Althoff, O. Stursberg, M. Buss.
"""
function expm_underapproximation(A::IntervalMatrix{T, Interval{T}}, t, p) where {T}
    n = checksquare(A)

    Y = zeros(n, n)
    LA = inf(A)
    Ail = LA * LA
    Z = zeros(n, n)
    RA = sup(A)
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

    B = IntervalMatrix(Matrix{Interval{T}}(undef, n , n))
    @inbounds for j in 1:n
        for i in 1:n
            minYZ = min(Y[i, j], Z[i, j])
            maxYZ = max(Y[i, j], Z[i, j])
            B[i, j] = Interval(minYZ, maxYZ)
        end
    end

    W = quadratic_expansion(A, t)
    res = W + B

    # add identity matrix implicitly
    for i in 1:n
        @inbounds res[i, i] += one(T)
    end
    return res
end

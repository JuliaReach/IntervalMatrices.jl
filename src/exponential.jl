using LinearAlgebra: checksquare

"""
    quadratic_expansion(A::IntervalMatrix, t)

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

    W = IntervalMatrix(similar(A.mat))

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
    quadratic_expansion(A::IntervalMatrix, α::Real, β::Real)

Compute the quadratic expansion of an interval matrix, ``αA + βA^2``, using
interval arithmetics.

### Input

- `A` -- interval matrix
- `α` -- linear coefficient
- `β` -- quadratic coefficient

### Output

An interval matrix that encloses ``B := αA + βA^2``.

### Algorithm

This a variation of the algorithm in [1, Section 6]. If ``A = (aᵢⱼ)`` and
``B = := αA + βA^2 = (bᵢⱼ)``, the idea is to compute each ``bᵢⱼ`` by factoring
out repeated expressions (thus the term *single-use expressions*).

First, let ``i ≡ j``. In this case,

```math
bᵢᵢ = β\\sum_\\{k, k ≠ i} a_{ik} a_{ki} + (α + βa_{ii}) a_{ii}.
```

Now consider ``i ≠ j``. Then,

```math
bᵢⱼ = β\\sum_\\{k, k ≠ i, k ≠ j} a_{ik} a_{kj} + (α + βa_{ii} + βa_{jj}) a_{ij}.
```

[1] Kosheleva, Kreinovich, Mayer, Nguyen. Computing the cube of an interval
matrix is NP-hard. SAC 2005.
"""
function quadratic_expansion(A::IntervalMatrix, α::Real, β::Real)
    B = similar(A.mat)
    n = checksquare(A)

    # case i == j
    @inbounds for j in 1:n
        B[j, j] = (α + β*A[j, j]) * A[j, j]
        for k in 1:n
            k == j && continue
            B[j, j] += β * (A[j, k] * A[k, j])
        end
    end

    # case i ≠ j
    @inbounds for j in 1:n
        for i in 1:n
            i == j && continue
            B[i, j] = A[i, j] * (α + β*A[j, j] + β*A[i, i])
            for k in 1:n
                (k == i || k == j) && continue
                B[i, j] += β * (A[i, k] * A[k, j])
            end
        end
    end
    return IntervalMatrix(B)
end

"""
    exp_overapproximation(M::IntervalMatrix{T, Interval{T}}, t, p) where {T}

Overapproximation of the exponential of an interval matrix.

### Input

- `A` -- interval matrix
- `t` -- non-negative time value
- `p` -- order of the approximation

### Algorithm

See Theorem 1 in *Reachability Analysis of Linear Systems with Uncertain
Parameters and Inputs* by M. Althoff, O. Stursberg, M. Buss.
"""
function exp_overapproximation(A::IntervalMatrix{T, Interval{T}}, t, p) where {T}
    n = checksquare(A)

    E = _exp_remainder(A, t, p; n=n)
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

# Implementation of Prop. 1 in Althoff, Matthias, Bruce H. Krogh, and Olaf Stursberg.
# "Analyzing reachability of linear dynamic systems with parametric uncertainties."
# Modeling, Design, and Simulation of Systems with Uncertainties. Springer, Berlin, Heidelberg, 2011. 69-94.
function _exp_remainder(A::IntervalMatrix{T}, t, p; n=checksquare(A)) where {T}
    C = max.(abs.(inf(A)), abs.(sup(A)))
    # compute Q = I + Ct + (Ct)^2/2! + ... + (Ct)^p/p!
    Q = Matrix(Diagonal(ones(T, n)))

    tⁱ = 1
    i! = 1
    Cⁱ = copy(C)
    for i in 1:p
        i! *= i
        tⁱ *= t
        Q = Q + Cⁱ * tⁱ/i!
        Cⁱ = Cⁱ * C
    end
    M = exp(C*t)
    Y  = M - Q
    Γ = IntervalMatrix(fill(zero(T)±one(T), (n, n)))
    E = Γ * Y
    return E
end

# Estimates the sum of the series in the matrix exponential. See Theorem 1
# in [1] Althoff, Matthias, Olaf Stursberg, and Martin Buss.
# Reachability analysis of nonlinear systems with uncertain parameters using conservative linearization.
# 2008 47th IEEE Conference on Decision and Control. IEEE, 2008.
function _exp_remainder_series(A::IntervalMatrix{T}, t, p; n=checksquare(A)) where {T}
    nA = opnorm(A, Inf)
    c = nA * t / (p + 2)
    @assert c < 1 "the remainder of the matrix exponential could not be " *
        "computed because a convergence condition is not satisfied: $c ≥ 1 " *
        "but it should be smaller than 1; try choosing a larger order"
    Γ = IntervalMatrix(fill(zero(T)±one(T), (n , n)))
    return Γ * ((nA*t)^(p+1) * (1/factorial(p + 1) * 1/(1-c)))
end

"""
    exp_underapproximation(M::IntervalMatrix{T, Interval{T}}, t, p) where {T}

Overapproximation of the exponential of an interval matrix.

### Input

- `A` -- interval matrix
- `t` -- non-negative time value
- `p` -- order of the approximation

### Algorithm

See Theorem 2 in *Reachability Analysis of Linear Systems with Uncertain
Parameters and Inputs* by M. Althoff, O. Stursberg, M. Buss.
"""
function exp_underapproximation(A::IntervalMatrix{T, Interval{T}}, t, p) where {T}
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

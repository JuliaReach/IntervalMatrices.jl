import Base: similar, split, ∈, ⊆, ∩
import Random: rand
import IntervalArithmetic: inf, sup, mid, diam

"""
    AbstractIntervalMatrix{IT} <: AbstractMatrix{IT}

Abstract supertype for interval matrix types.
"""
abstract type AbstractIntervalMatrix{IT} <: AbstractMatrix{IT} end

"""
    IntervalMatrix{T, IT<:Interval{T}, MT<:AbstractMatrix{IT}} <: AbstractIntervalMatrix{IT}

An interval matrix i.e. a matrix whose coefficients are intervals. This type is
parametrized in the number field, the interval type, and the matrix type.

### Fields

- `mat` -- matrix whose entries are intervals

### Examples

```jldoctest
julia> A = IntervalMatrix([-0.9±0.1 0±0; 0±0 -0.9±0.1])
2×2 IntervalMatrix{Float64,Interval{Float64},Array{Interval{Float64},2}}:
 [-1, -0.799999]   [0, 0]
  [0, 0]          [-1, -0.799999]
```

An interval matrix proportional to the identity matrix can be built using the
`UniformScaling` operator from the standard library `LinearAlgebra`. For example,

```jldoctest interval_uniform_scaling
julia> using LinearAlgebra

julia> IntervalMatrix(Interval(1)*I, 2)
2×2 IntervalMatrix{Float64,Interval{Float64},Array{Interval{Float64},2}}:
 [1, 1]  [0, 0]
 [0, 0]  [1, 1]
```
The number of columns can be specified as a third argument, creating a rectangular
``m × n`` matrix such that only the entries in the main diagonal,
``(1, 1), (2, 2), …,  (k, k)`` are specified, where ``k = \\min(m, n)``:

```jldoctest interval_uniform_scaling
julia> IntervalMatrix(Interval(-1, 1)*I, 2, 3)
2×3 IntervalMatrix{Float64,Interval{Float64},Array{Interval{Float64},2}}:
 [-1, 1]   [0, 0]  [0, 0]
  [0, 0]  [-1, 1]  [0, 0]

julia> IntervalMatrix(Interval(-1, 1)*I, 3, 2)
3×2 IntervalMatrix{Float64,Interval{Float64},Array{Interval{Float64},2}}:
 [-1, 1]   [0, 0]
  [0, 0]  [-1, 1]
  [0, 0]   [0, 0]
```

An uninitialized interval matrix can be constructed using `undef`:

```jldoctest undef_test
julia> m = IntervalMatrix{Float64}(undef, 2, 2);

julia> typeof(m)
IntervalMatrix{Float64,Interval{Float64},Array{Interval{Float64},2}}
```
Note that this constructor implicitly uses a dense matrix, `Matrix{Float64}`,
as the matrix (`mat`) field in the new interval matrix.
"""
struct IntervalMatrix{T, IT<:Interval{T}, MT<:AbstractMatrix{IT}} <: AbstractIntervalMatrix{IT}
    mat::MT
end

import Base:size, IndexStyle, getindex, setindex!, copy

IndexStyle(::Type{<:IntervalMatrix}) = IndexLinear()
size(M::IntervalMatrix) = size(M.mat)
getindex(M::IntervalMatrix, i::Int) = getindex(M.mat, i)
setindex!(M::IntervalMatrix, X, inds...) = setindex!(M.mat, X, inds...)
copy(M::IntervalMatrix) = IntervalMatrix(copy(M.mat))

# constructor from uniform scaling
function IntervalMatrix(αI::UniformScaling{Interval{T}}, m::Integer, n::Integer=m) where {T}
    return IntervalMatrix(Matrix(αI, m, n))
end

# undef initializer, eg. IntervalMatrix{Float64}(undef, 2, 2)
function IntervalMatrix{T}(u::UndefInitializer, m::Integer, n::Integer=m) where {T}
    mat = Matrix{Interval{T}}(undef, m, n)
    return IntervalMatrix(mat)
end

# similar initializer
function similar(M::IntervalMatrix)
    return IntervalMatrix(similar(M.mat))
end

# constructor from a scalar matrix
function IntervalMatrix(A::AbstractMatrix{T}) where {T<:Number}
    return IntervalMatrix(map(Interval, A))
end

"""
    IntervalMatrix(C::MT, S::MT) where {T, MT<:AbstractMatrix{T}}

Return an interval matrix such that the center and radius of the intervals
is given by the matrices `C` and `S` respectively.

### Input

- `C` -- center matrix
- `S` -- radii matrix

### Output

An interval matrix `M` such that `M[i, j]` corresponds to the interval whose center
is `C[i, j]` and whose radius is `S[i, j]`, for each `i` and `j`. That is,
``M = C + [-S, S]``.

### Notes

The radii matrix should be nonnegative, i.e. `S[i, j] ≥ 0` for each `i` and `j`.

### Examples

```jldoctest
julia> IntervalMatrix([1 2; 3 4], [1 2; 4 5])
2×2 IntervalMatrix{Float64,Interval{Float64},Array{Interval{Float64},2}}:
  [0, 2]   [0, 4]
 [-1, 7]  [-1, 9]
 ```
"""
function IntervalMatrix(C::MT, S::MT) where {T, MT<:AbstractMatrix{T}}
    size(C) == size(S) || throw(ArgumentError("the sizes of the center matrix and the " *
                                "radii matrix should match, but they are $(size(C)) " *
                                "and $(size(S)) respectively"))
    m, n = size(C)
    M = IntervalMatrix{T}(undef, m, n)

    @inbounds for j in 1:n
        for i in 1:m
            M[i, j] = C[i, j] ± S[i, j]
        end
    end

    return M
end

"""
    opnorm(A::IntervalMatrix, p::Real=Inf)

The matrix norm of an interval matrix.

### Input

- `A` -- interval matrix
- `p` -- (optional, default: `Inf`) the class of `p`-norm

### Notes

The matrix ``p``-norm of an interval matrix ``A`` is defined as

```math
    ‖A‖_p := ‖\\max(|\\text{inf}(A)|, |\\text{sup}(A)|)‖_p
```

where ``\\max`` and ``|·|`` are taken elementwise.
"""
function LinearAlgebra.opnorm(A::IntervalMatrix, p::Real=Inf)
    if p == Inf
        return _opnorm_inf(A)
    elseif p == 1
        return _opnorm_1(A)
    else
        error("the interval matrix norm for this value of p=$p is not implemented")
    end
end

# The operator norm in the infinity norm corresponds to the
# maximum absolute value row-sum.
function _opnorm_inf(A::IntervalMatrix{T}) where {T}
    m, n = size(A)
    res = zero(T)
    @inbounds @simd for i in 1:m
        acc = zero(T)
        for j in 1:n
            x = A[i, j]
            acc += max(abs(inf(x)), abs(sup(x)))
        end
        if acc > res
            res = acc
        end
    end
    return res
end

# The operator norm in the 1-norm corresponds to the
# maximum absolute value column-sum.
function _opnorm_1(A::IntervalMatrix{T}) where {T}
    m, n = size(A)
    res = zero(T)
    @inbounds @simd for j in 1:n
        acc = zero(T)
        for i in 1:m
            x = A[i, j]
            acc += max(abs(inf(x)), abs(sup(x)))
        end
        if acc > res
            res = acc
        end
    end
    return res
end

"""
    inf(A::IntervalMatrix{T}) where {T}

Return the infimum of an interval matrix `A`, which corresponds to taking the
element-wise infimum of `A`.

### Input

- `A` -- interval matrix

### Output

A scalar matrix whose coefficients are the infima of each element in `A`.
"""
function inf(A::IntervalMatrix{T}) where {T}
    return map(inf, A)
end

"""
    sup(A::IntervalMatrix{T}) where {T}

Return the supremum of an interval matrix `A`, which corresponds to taking the
element-wise supremum of `A`.

### Input

- `A` -- interval matrix

### Output

A scalar matrix whose coefficients are the suprema of each element in `A`.
"""
function sup(A::IntervalMatrix{T}) where {T}
    return map(sup, A)
end

"""
    mid(A::IntervalMatrix{T}) where {T}

Return the midpoint of an interval matrix `A`, which corresponds to taking the
element-wise midpoint of `A`.

### Input

- `A` -- interval matrix

### Output

A scalar matrix whose coefficients are the midpoints of each element in `A`.
"""
function mid(A::IntervalMatrix{T}) where {T}
    return map(mid, A)
end

"""
    split(A::IntervalMatrix{T}) where {T}

Split an interval matrix ``A`` into two scalar matrices ``C`` and ``S``
such that ``A = C + [-S, S]``.

### Input

- `A` -- interval matrix

### Output

A pair `(C, S)` such that the entries of `C` are the central points and the
entries of `S` are the (nonnegative) radii of the intervals in `A`.
"""
function split(A::IntervalMatrix{T}) where {T}
    m, n = size(A)
    C = Matrix{T}(undef, m, n)
    S = Matrix{T}(undef, m, n)

    @inbounds for j in 1:n
        for i in 1:m
            itv = A[i, j]
            radius = (sup(itv) - inf(itv)) / T(2)
            C[i, j] = inf(itv) + radius
            S[i, j] = radius
        end
    end

    return C, S
end

"""
    ∈(M::AbstractMatrix, A::AbstractIntervalMatrix)

Check whether a concrete matrix is an instance of an interval matrix.

### Input

- `M` -- concrete matrix
- `A` -- interval matrix

### Output

`true` iff `M` is an instance of `A`

### Algorithm

We check for each entry in `M` whether it belongs to the corresponding interval
in `A`.
"""
function ∈(M::AbstractMatrix, A::AbstractIntervalMatrix)
    @assert size(M) == size(A) "incompatible matrix sizes (M: $(size(M)), A: " *
                               "$(size(A)))"

    m, n = size(A)
    @inbounds for j in 1:n
        for i in 1:m
            if M[i, j] ∉ A[i, j]
                return false
            end
        end
    end
    return true
end

"""
    ⊆(A::AbstractIntervalMatrix, B::AbstractIntervalMatrix)

Check whether an interval matrix is contained in another interval matrix.

### Input

- `A` -- interval matrix
- `B` -- interval matrix

### Output

`true` iff `A[i, j] ⊆ B[i, j]` for all `i, j`.
"""
function ⊆(A::AbstractIntervalMatrix, B::AbstractIntervalMatrix)
    @assert size(A) == size(B) "incompatible matrix sizes $(size(A)) and " *
                               "$(size(B))"

    m, n = size(A)
    @inbounds for j in 1:n, i in 1:m
        if !(A[i, j] ⊆ B[i, j])
            return false
        end
    end
    return true
end

# random interval
function rand(::Type{Interval}; N::Type{<:Real}=Float64,
              rng::AbstractRNG=GLOBAL_RNG)::Interval{N}
    x, y = randn(rng, N), randn(rng, N)
    return x < y ? Interval(x, y) : Interval(y, x)
end

# create a random interval for the given numeric type and random-number generator
@inline function _rand_interval(; N=Float64, rng::AbstractRNG=GLOBAL_RNG)
    x, y = randn(rng, N), randn(rng, N)
    return x < y ? Interval(x, y) : Interval(y, x)
end

"""
    rand(::Type{IntervalMatrix}, m::Int=2, [n]::Int=m;
         N=Float64, rng::AbstractRNG=GLOBAL_RNG)

Return a random interval matrix of the given size and numeric type.

### Input

- `IntervalMatrix` -- type, used for dispatch
- `m`              -- (optional, default: `2`) number of rows
- `n`              -- (optional, default: `m`) number of columns
- `rng`            -- (optional, default: `GLOBAL_RNG`) random-number generator

### Output

An interval matrix of size ``m × n`` whose coefficients are normally-distributed
intervals of type `N` with mean `0` and standard deviation `1`.

### Notes

If this function is called with only one argument, it creates a square matrix,
because the number of columns defaults to the number of rows.
"""
function rand(::Type{IntervalMatrix}, m::Int=2, n::Int=m;
              N=Float64, rng::AbstractRNG=GLOBAL_RNG)
    B = Matrix{Interval{N}}(undef, m, n)
    for j in 1:n
        for i in 1:m
            @inbounds B[i, j] = _rand_interval(N=N, rng=rng)
        end
    end
    return IntervalMatrix(B)
end

"""
    sample(A::IntervalMatrix{T}; rng::AbstractRNG=GLOBAL_RNG) where {T}

Return a sample of the given random interval matrix.

### Input

- `A`   -- interval matrix
- `m`   -- (optional, default: `2`) number of rows
- `n`   -- (optional, default: `2`) number of columns
- `rng` -- (optional, default: `GLOBAL_RNG`) random-number generator

### Output

An interval matrix of size ``m × n`` whose coefficients are normally-distributed
intervals of type `N` with mean `0` and standard deviation `1`.
"""
function sample(A::IntervalMatrix{T}; rng::AbstractRNG=GLOBAL_RNG) where {T}
    m, n = size(A)
    B = Matrix{T}(undef, m, n)
    @inbounds for j in 1:n
        for i in 1:m
            itv = A[i, j]
            B[i, j] = (sup(itv) - inf(itv)) * rand(rng) + inf(itv)
        end
    end
    return B
end

"""
    diam(A::IntervalMatrix{T}) where {T}

Return a matrix whose entries describe the diameters of the intervals.

### Input

- `A` -- interval matrix

### Output

A matrix `B` of the same shape as `A` such that `B[i, j] == diam(A[i, j])` for
each `i` and `j`.
"""
function diam(A::IntervalMatrix{T}) where {T}
    return map(diam, A)
end

"""
    scale(A::IntervalMatrix{T}, α::T) where {T}

Return a new interval matrix whose entries are scaled by the given factor.

### Input

- `A` -- interval matrix
- `α` -- scaling factor

### Output

A new matrix `B` of the same shape as `A` such that `B[i, j] = α*A[i, j]` for
each `i` and `j`.

### Notes

See `scale!` for the in-place version of this function.
"""
function scale(A::IntervalMatrix{T}, α::T) where {T}
    return scale!(copy(A), α)
end

"""
    scale(A::IntervalMatrix{T}, α::T) where {T}

Modifies the given interval matrix, scaling its entries by the given factor.

### Input

- `A` -- interval matrix
- `α` -- scaling factor

### Output

The matrix `A` such that for each `i` and `j`, the new value of `A[i, j]` is
`α*A[i, j]`.

### Notes

This is the in-place version of `scale`.
"""
function scale!(A::IntervalMatrix{T}, α::T) where {T}
    return map!(x -> α * x, A, A)
end

"""
    ∩(A::IntervalMatrix, B::IntervalMatrix)

Intersect two interval matrices.

### Input

- `A` -- interval matrix
- `B` -- interval matrix (of the same shape as `A`)

### Output

A new matrix `C` of the same shape as `A` such that
`C[i, j] = A[i, j] ∩ B[i, j]` for each `i` and `j`.
"""
function ∩(A::IntervalMatrix, B::IntervalMatrix)
    m, n = size(A)
    @assert size(A) == size(B) "incompatible matrix sizes (A: $(size(A)), B: " *
                               "$(size(B)))"

    C = similar(A)
    @inbounds for j in 1:n, i in 1:m
        C[i, j] = A[i, j] ∩ B[i, j]
    end
    return C
end

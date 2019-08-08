import Base: split,
             ∈
import Random: rand

"""
    AbstractIntervalMatrix{IT} <: AbstractMatrix{IT}

Abstract supertype for interval matrix types.
"""
abstract type AbstractIntervalMatrix{IT} <: AbstractMatrix{IT} end

"""
    IntervalMatrix{T, IT<:Interval{T}, MT<:AbstractMatrix{IT}} <: AbstractIntervalMatrix{IT}

An interval matrix i.e. a matrix whose coefficients are intervals. This type is
parametrized in the number field, the type of interval.

### Fields

- `mat` -- matrix whose entries are intervals

### Examples

```jldoctest
julia> A = IntervalMatrix([-0.9±0.1 0±0; 0±0 -0.9±0.1])
2×2 IntervalMatrix{Float64,Interval{Float64},Array{Interval{Float64},2}}:
 [-1, -0.799999]   [0, 0]
  [0, 0]          [-1, -0.799999]
```
"""
struct IntervalMatrix{T, IT<:Interval{T}, MT<:AbstractMatrix{IT}} <: AbstractIntervalMatrix{IT}
    mat::MT

    # default constructor with domain constraint for radius
    function IntervalMatrix{T, IT, MT}(M::MT) where {T, IT<:Interval{T}, MT<:AbstractMatrix{IT}}
        return new{T, IT, MT}(M)
    end
end

import Base:size, IndexStyle, getindex, setindex!

IndexStyle(::Type{<:IntervalMatrix}) = IndexLinear()
size(M::IntervalMatrix) = size(M.mat)
getindex(M::IntervalMatrix, i::Int) = getindex(M.mat, i)
setindex!(M::IntervalMatrix, X, inds...) = setindex!(M.mat, X, inds...)

# convenience constructor without type parameter
function IntervalMatrix(M::MT) where {T, IT<:Interval{T}, MT<:AbstractMatrix{IT}}
    return IntervalMatrix{T, IT, MT}(M)
end

# convenience constructor for undef initializer
function IntervalMatrix{T, IT, MT}(::UndefInitializer, m::Int, n::Int) where {T, IT<:Interval{T}, MT<:AbstractMatrix{IT}}
    return IntervalMatrix(MT(undef, m, n))
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
    ‖A‖_p := ‖\\max(|\\text{left}(A)|, |\\text{right}(A)|)‖_p
```

where ``\\max`` and ``|·|`` are taken elementwise.
"""
function LinearAlgebra.opnorm(A::IntervalMatrix, p::Real=Inf)
    if p == Inf || p == 1
        return LinearAlgebra.opnorm(max.(abs.(left(A)), abs.(right(A))), p)
    else
        error("the interval matrix norm for this value of p=$p is not implemented")
    end
end

"""
    left(A::IntervalMatrix{T}) where {T}

Return the left part of this interval matrix, which corresponds to taking the
element-wise infimum of `A`.

### Input

- `A` -- interval matrix

### Output

A scalar matrix whose coefficients are the infima of each element in `A`.
"""
function left(A::IntervalMatrix{T}) where {T}
    m, n = size(A)
    L = Matrix{T}(undef, m, n)
    for j in 1:m
        for i in 1:n
            @inbounds L[i, j] = inf(A[i, j])
        end
    end
    return L
end

"""
    right(A::IntervalMatrix{T}) where {T}

Return the right part of this interval matrix, which corresponds to taking the
element-wise supremum of `A`.

### Input

- `A` -- interval matrix

### Output

A scalar matrix whose coefficients are the suprema of each element in `A`.
"""
function right(A::IntervalMatrix{T}) where {T}
    m, n = size(A)
    R = Matrix{T}(undef, m, n)
    for j in 1:m
        for i in 1:n
            @inbounds R[i, j] = sup(A[i, j])
        end
    end
    return R
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

# random interval
function rand(::Type{Interval}; N::Type{<:Real}=Float64,
              rng::AbstractRNG=GLOBAL_RNG)::Interval{N}
    x, y = randn(rng, N), randn(rng, N)
    return x < y ? Interval(x, y) : Interval(y, x)
end

"""
    rand(::Type{IntervalMatrix}, m::Int=2, n::Int=2;
         N=Float64, rng::AbstractRNG=GLOBAL_RNG)

Return a random interval matrix of the given size and numeric type.

### Input

- `IntervalMatrix` -- type, used for dispatch
- `m`              -- (optional, default: `2`) number of rows
- `n`              -- (optional, default: `2`) number of columns
- `rng`            -- (optional, default: `GLOBAL_RNG`) random-number generator

### Output

An interval matrix of size ``m × n`` whose coefficients are normally-distributed
intervals of type `N` with mean `0` and standard deviation `1`.
"""
function rand(::Type{IntervalMatrix}, m::Int=2, n::Int=2;
              N=Float64, rng::AbstractRNG=GLOBAL_RNG)
    B = Matrix{Interval{N}}(undef, m, n)
    for j in 1:n
        for i in 1:m
            @inbounds B[i, j] = rand(Interval, N=N)
        end
    end
    return IntervalMatrix(B)
end

"""
    sample(A::IntervalMatrix{T}; rng::AbstractRNG=GLOBAL_RNG) where {T}

Return a sample of the given random interval matrix.

### Input

- `A` -- interval matrix
- `m`              -- (optional, default: `2`) number of rows
- `n`              -- (optional, default: `2`) number of columns
- `rng`            -- (optional, default: `GLOBAL_RNG`) random-number generator

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

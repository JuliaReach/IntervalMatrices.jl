"""
    AbstractIntervalMatrix{IT} <: AbstractMatrix{IT}

Abstract supertype for interval matrix types.
"""
abstract type AbstractIntervalMatrix{IT} <: AbstractMatrix{IT} end

"""
    IntervalMatrix{T, IT<:AbstractInterval{T}, M<:AbstractMatrix{IT}} <:
        AbstractIntervalMatrix{IT}

Interval matrix, this type is parametrized in the number field and the type of
interval set.

### Fields

- `mat` -- matrix whose entries are intervals

### Examples

```jldoctest
julia> A = IntervalMatrix([-0.9±0.1 0±0; 0±0 -0.9±0.1])
2×2 IntervalMatrix{Float64,Interval{:closed,:closed,Float64},Array{Interval{:closed,:closed,Float64},2}}:
 -1.0..-0.8  0.0..0.0
 0.0..0.0    -1.0..-0.8
```
"""
struct IntervalMatrix{T, IT<:AbstractInterval{T}, M<:AbstractMatrix{IT}} <:
        AbstractIntervalMatrix{IT}
    mat::M
end

import Base:size, IndexStyle, getindex, setindex!, +, *

IndexStyle(::Type{<:IntervalMatrix}) = IndexLinear()
size(M::IntervalMatrix) = size(M.mat)
getindex(M::IntervalMatrix, i::Int) = getindex(M.mat, i)
setindex!(M::IntervalMatrix, X, inds...) = setindex!(M.mat, X, inds...)

+(M1::IntervalMatrix, M2::IntervalMatrix) = IntervalMatrix(M1.mat + M2.mat)
*(M1::IntervalMatrix, M2::IntervalMatrix) = IntervalMatrix(M1.mat * M2.mat)

"""
    opnorm(A::IntervalMatrix, p::Real=Inf)

The matrix norm of an interval matrix.

### Input

- `A` -- interval matrix
- `p` -- (optional, default: `Inf`): the class of `p`-norm

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
    left(A::IntervalMatrix)

The left part of this interval matrix, which corresponds to the `left` of
each interval element in the matrix.

### Input

- `A` -- interval matrix
"""
function left(A::IntervalMatrix)
    n = size(A, 1)
    return hcat([[A[i, j].left for j in 1:n] for i in 1:n]...)'
end

"""
    right(A::IntervalMatrix)

The right part of this interval matrix, which corresponds to the `right` of
each interval element in the matrix.

### Input

- `A` -- interval matrix
"""
function right(A::IntervalMatrix)
    n = size(A, 1)
    return hcat([[A[i, j].right for j in 1:n] for i in 1:n]...)'
end

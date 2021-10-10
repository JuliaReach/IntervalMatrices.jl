"""
    IntervalMatrixPencil{T, IT, MT0<:AbstractMatrix{T}, MT1<:AbstractMatrix{T}} <: AbstractIntervalMatrix{IT}

Interval matrix representing the matrix

```math
A₀ + λA₁,
```
where ``A₀`` and ``A₁`` are real (or complex) matrices, and ``λ`` is an interval.

### Fields

- `A0` -- matrix
- `A1` -- matrix
- `λ`  -- interval

### Examples

The matrix pencil ``I + [1 1; -1 1] * (0 .. 1)`` is:

```jldoctest
julia> using LinearAlgebra

julia> P = IntervalMatrixPencil(Matrix(1.0I, 2, 2), [1 1; -1 1.], 0 .. 1);

julia> P
2×2 IntervalMatrixPencil{Float64, Interval{Float64}, Matrix{Float64}, Matrix{Float64}}:
  [1, 2]  [0, 1]
 [-1, 0]  [1, 2]
```
"""
struct IntervalMatrixPencil{T, IT, MT0<:AbstractMatrix{T}, MT1<:AbstractMatrix{T}} <: AbstractIntervalMatrix{IT}
    A0::MT0
    A1::MT1
    λ::IT

    # inner constructor with dimension check
    function IntervalMatrixPencil(A0::MT0, A1::MT1, λ::IT) where {T, IT, MT0<:AbstractMatrix{T}, MT1<:AbstractMatrix{T}}
        @assert checksquare(A0) == checksquare(A1) "the size of `A0` and `A1` should match, got $(size(A0)) and $(size(A1)) respectively"
        return new{T, IT, MT0, MT1}(A0, A1, λ)
    end
end

IndexStyle(::Type{<:IntervalMatrixPencil}) = IndexLinear()
size(M::IntervalMatrixPencil) = size(M.A0)
getindex(M::IntervalMatrixPencil, i::Int) = getindex(M.A0, i) + M.λ * getindex(M.A1, i)
function setindex!(M::IntervalMatrixPencil{T}, X::T, inds...) where {T}
    setindex!(M.A0, X, inds...)
    if !iszero(getindex(M.A1, inds...)) || !iszero(M.λ)
        setindex!(M.A1, zero(T), inds...)
    end
end
copy(M::IntervalMatrixPencil) = IntervalMatrixPencil(copy(M.A0), copy(M.A1), M.λ)

"""
    AffineIntervalMatrix{T, IT, MT0<:AbstractMatrix{T}, MT<:AbstractMatrix{T}, MTA<:AbstractVector{MT}} <: AbstractIntervalMatrix{IT}

Interval matrix representing the matrix

```math
A₀ + λ₁A₁ + λ₂A₂ + … + λₖAₖ,
```
where ``A₀`` and ``A₁, …, Aₖ`` are real (or complex) matrices, and ``λ₁, …, λₖ``
are intervals.

### Fields

- `A0` -- matrix
- `A`  -- vector of matrices
- `λ`  -- vector of intervals

### Notes

This type is the general case of the [`IntervalMatrixPencil`](@ref), which only
contains one matrix proportional to an interval.

### Examples

The matrix pencil ``I + [1 1; -1 1] * (0 .. 1) + [0 1; 1 0] * (2 .. 3)`` is:

```jldoctest
julia> using LinearAlgebra

julia> A0 = Matrix(1.0I, 2, 2);

julia> A1 = [1 1; -1 1.]; A2 = [0 1; 1 0];

julia> λ1 = 0 .. 1; λ2 = 2 .. 3;

julia> P = AffineIntervalMatrix(A0, [A1, A2], [λ1, λ1]);

julia> P[1, 1]
[1, 3]

julia> P[1, 2]
[0, 2]
```
"""
struct AffineIntervalMatrix{T, IT, MT0<:AbstractMatrix{T}, MT<:AbstractMatrix{T}, MTA<:AbstractVector{MT}, VIT<:AbstractVector{IT}} <: AbstractIntervalMatrix{IT}
    A0::MT0
    A::MTA
    λ::VIT
end

IndexStyle(::Type{<:AffineIntervalMatrix}) = IndexLinear()
size(M::AffineIntervalMatrix) = size(M.A0)
getindex(M::AffineIntervalMatrix, i::Int) = getindex(M.A0, i) + sum(M.λ[k] * getindex(M.A[k], i) for k in eachindex(M.λ))
function setindex!(M::AffineIntervalMatrix{T}, X::T, inds...) where {T}
    setindex!(M.A0, X, inds...)
    @inbounds for k in 1:length(M.A)
        setindex!(M.A[k], zero(T), inds...)
    end
end
copy(M::AffineIntervalMatrix) = AffineIntervalMatrix(copy(M.A0), copy(M.A), M.λ)

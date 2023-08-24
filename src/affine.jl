"""
    AffineIntervalMatrix1{T, IT, MT0<:AbstractMatrix{T}, MT1<:AbstractMatrix{T}} <: AbstractIntervalMatrix{IT}

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

The matrix ``I + [1 1; -1 1] * (0 .. 1)`` is:

```jldoctest
julia> using LinearAlgebra

julia> P = AffineIntervalMatrix1(Matrix(1.0I, 2, 2), [1 1; -1 1.], 0 .. 1);

julia> P
2×2 AffineIntervalMatrix1{Float64, Interval{Float64}, Matrix{Float64}, Matrix{Float64}}:
  [1.0, 2.0]  [0.0, 1.0]
 [-1.0, 0.0]  [1.0, 2.0]
```
"""
struct AffineIntervalMatrix1{T,IT,MT0<:AbstractMatrix{T},MT1<:AbstractMatrix{T}} <:
       AbstractIntervalMatrix{IT}
    A0::MT0
    A1::MT1
    λ::IT

    # inner constructor with dimension check
    function AffineIntervalMatrix1(A0::MT0, A1::MT1,
                                   λ::IT) where {T,IT,MT0<:AbstractMatrix{T},MT1<:AbstractMatrix{T}}
        @assert checksquare(A0) == checksquare(A1) "the size of `A0` and `A1` should match, got $(size(A0)) and $(size(A1)) respectively"
        return new{T,IT,MT0,MT1}(A0, A1, λ)
    end
end

IndexStyle(::Type{<:AffineIntervalMatrix1}) = IndexLinear()
size(M::AffineIntervalMatrix1) = size(M.A0)
getindex(M::AffineIntervalMatrix1, i::Int) = getindex(M.A0, i) + M.λ * getindex(M.A1, i)
function setindex!(M::AffineIntervalMatrix1{T}, X::T, inds...) where {T}
    setindex!(M.A0, X, inds...)
    return setindex!(M.A1, zero(T), inds...)
end
copy(M::AffineIntervalMatrix1) = AffineIntervalMatrix1(copy(M.A0), copy(M.A1), M.λ)

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

This type is the general case of the [`AffineIntervalMatrix1`](@ref), which only
contains one matrix proportional to an interval.

### Examples

The affine matrix ``I + [1 1; -1 1] * (0 .. 1) + [0 1; 1 0] * (2 .. 3)`` is:

```jldoctest
julia> using LinearAlgebra

julia> A0 = Matrix(1.0I, 2, 2);

julia> A1 = [1 1; -1 1.]; A2 = [0 1; 1 0];

julia> λ1 = 0 .. 1; λ2 = 2 .. 3;

julia> P = AffineIntervalMatrix(A0, [A1, A2], [λ1, λ1]);

julia> P[1, 1]
[1.0, 3.0]

julia> P[1, 2]
[0.0, 2.0]
```
"""
struct AffineIntervalMatrix{T,IT,MT0<:AbstractMatrix{T},MT<:AbstractMatrix{T},
                            MTA<:AbstractVector{MT},VIT<:AbstractVector{IT}} <:
       AbstractIntervalMatrix{IT}
    A0::MT0
    A::MTA
    λ::VIT

    # inner constructor with dimension check
    function AffineIntervalMatrix(A0::MT0, A::MTA,
                                  λ::VIT) where {T,IT,MT0<:AbstractMatrix{T},MT<:AbstractMatrix{T},
                                                 MTA<:AbstractVector{MT},VIT<:AbstractVector{IT}}
        k = length(A)
        @assert k == length(λ) "expected `A` and `λ` to have the same length, got $(length(A)) and $(length(λ)) respectively"
        n = checksquare(A0)
        @inbounds for i in 1:k
            @assert n == checksquare(A[i]) "each matrix should have the same size $n × $n"
        end
        return new{T,IT,MT0,MT,MTA,VIT}(A0, A, λ)
    end
end

IndexStyle(::Type{<:AffineIntervalMatrix}) = IndexLinear()
size(M::AffineIntervalMatrix) = size(M.A0)
function getindex(M::AffineIntervalMatrix, i::Int)
    return getindex(M.A0, i) + sum(M.λ[k] * getindex(M.A[k], i) for k in eachindex(M.λ))
end
function setindex!(M::AffineIntervalMatrix{T}, X::T, inds...) where {T}
    setindex!(M.A0, X, inds...)
    @inbounds for k in 1:length(M.A)
        setindex!(M.A[k], zero(T), inds...)
    end
end
copy(M::AffineIntervalMatrix) = AffineIntervalMatrix(copy(M.A0), copy(M.A), M.λ)

"""
    AbstractIntervalMatrix{IT} <: AbstractMatrix{IT}

Abstract supertype for interval matrix types.
"""
abstract type AbstractIntervalMatrix{IT} <: AbstractMatrix{IT} end

"""
    IntervalMatrix{T, IT, MT<:AbstractMatrix{IT}} <: AbstractIntervalMatrix{IT}

An interval matrix i.e. a matrix whose coefficients are intervals. This type is
parameterized in the number field, the interval type, and the matrix type.

### Fields

- `mat` -- matrix whose entries are intervals

### Examples

```jldoctest
julia> A = IntervalMatrix([interval(-1, -0.8) interval(0); interval(0) interval(-1, -0.8)])
2×2 IntervalMatrix{Float64, Interval{Float64}, Matrix{Interval{Float64}}}:
 [-1.0, -0.8]   [0.0, 0.0]
  [0.0, 0.0]   [-1.0, -0.8]
```

An interval matrix proportional to the identity matrix can be built using the
`UniformScaling` operator from the standard library `LinearAlgebra`. For example,

```jldoctest interval_uniform_scaling
julia> using LinearAlgebra

julia> IntervalMatrix(interval(1)*I, 2)
2×2 IntervalMatrix{Float64, Interval{Float64}, Matrix{Interval{Float64}}}:
 [1.0, 1.0]  [0.0, 0.0]
 [0.0, 0.0]  [1.0, 1.0]
```
The number of columns can be specified as a third argument, creating a rectangular
``m × n`` matrix such that only the entries in the main diagonal,
``(1, 1), (2, 2), …,  (k, k)`` are specified, where ``k = \\min(m, n)``:

```jldoctest interval_uniform_scaling
julia> IntervalMatrix(interval(-1, 1)*I, 2, 3)
2×3 IntervalMatrix{Float64, Interval{Float64}, Matrix{Interval{Float64}}}:
 [-1.0, 1.0]   [0.0, 0.0]  [0.0, 0.0]
  [0.0, 0.0]  [-1.0, 1.0]  [0.0, 0.0]

julia> IntervalMatrix(interval(-1, 1)*I, 3, 2)
3×2 IntervalMatrix{Float64, Interval{Float64}, Matrix{Interval{Float64}}}:
 [-1.0, 1.0]   [0.0, 0.0]
  [0.0, 0.0]  [-1.0, 1.0]
  [0.0, 0.0]   [0.0, 0.0]
```

An uninitialized interval matrix can be constructed using `undef`:

```jldoctest undef_test
julia> m = IntervalMatrix{Float64}(undef, 2, 2);

julia> typeof(m)
IntervalMatrix{Float64, Interval{Float64}, Matrix{Interval{Float64}}}
```
Note that this constructor implicitly uses a dense matrix, `Matrix{Float64}`,
as the matrix (`mat`) field in the new interval matrix.
"""
struct IntervalMatrix{T,IT,MT<:AbstractMatrix{IT}} <: AbstractIntervalMatrix{IT}
    mat::MT
end

IndexStyle(::Type{<:IntervalMatrix}) = IndexLinear()
size(M::IntervalMatrix) = size(M.mat)
getindex(M::IntervalMatrix, i::Int) = getindex(M.mat, i)
setindex!(M::IntervalMatrix, X, inds...) = setindex!(M.mat, X, inds...)
copy(M::IntervalMatrix) = IntervalMatrix(copy(M.mat))

function IntervalMatrix(A::MT) where {T,IT<:Interval{T},MT<:AbstractMatrix{IT}}
    return IntervalMatrix{T,IT,MT}(A)
end

function IntervalMatrix(A::MT) where {T,IT<:Complex{Interval{T}},MT<:AbstractMatrix{IT}}
    return IntervalMatrix{T,IT,MT}(A)
end

# constructor from uniform scaling
function IntervalMatrix(αI::UniformScaling{Interval{T}}, m::Integer, n::Integer=m) where {T}
    return IntervalMatrix(Matrix(αI, m, n))
end

# undef initializer, eg. IntervalMatrix{Float64}(undef, 2, 2)
function IntervalMatrix{T}(u::UndefInitializer, m::Integer, n::Integer=m) where {T}
    mat = Matrix{Interval{T}}(u, m, n)
    return IntervalMatrix(mat)
end

# similar initializer
function similar(M::IntervalMatrix)
    return IntervalMatrix(similar(M.mat))
end

# constructor from a scalar matrix
function IntervalMatrix(A::AbstractMatrix{T}) where {T<:Number}
    return IntervalMatrix(map(interval, A))
end

"""
    IntervalMatrix(A::MT, B::MT) where {T, MT<:AbstractMatrix{T}}

Return an interval matrix such that the lower and upper bounds of the intervals
are given by the matrices `A` and `B` respectively.

### Input

- `A` -- lower bound matrix
- `B` -- upper bound matrix

### Output

An interval matrix `M` such that `M[i, j]` corresponds to the interval whose lower bound
is `A[i, j]` and whose upper bound is `B[i, j]`, for each `i` and `j`. That is,
``M_{ij} = [A_{ij}, B_{ij}]``.

### Notes

The upper bound should be bigger or equal than the lower bound,
i.e. `B[i, j] ≥ A[i, j]` for each `i` and `j`.

### Examples

```jldoctest
julia> IntervalMatrix([1 2; 3 4], [1 2; 4 5])
2×2 IntervalMatrix{Float64, Interval{Float64}, Matrix{Interval{Float64}}}:
 [1.0, 1.0]  [2.0, 2.0]
 [3.0, 4.0]  [4.0, 5.0]
```
"""
function IntervalMatrix(A::MT, B::MT) where {T,MT<:AbstractMatrix{T}}
    size(A) == size(B) || throw(ArgumentError("the sizes of the lower and upper bound " *
                                              "matrices should match, but they are $(size(A)) " *
                                              "and $(size(B)) respectively"))

    return IntervalMatrix(map((x, y) -> interval(x, y), A, B))
end

"""
    ±(C::MT, S::MT) where {T, MT<:AbstractMatrix{T}}

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
julia> [1 2; 3 4] ± [1 2; 4 5]
2×2 IntervalMatrix{Float64, Interval{Float64}, Matrix{Interval{Float64}}}:
  [0.0, 2.0]   [0.0, 4.0]
 [-1.0, 7.0]  [-1.0, 9.0]
```
"""
function ±(C::MT, S::MT) where {T,MT<:AbstractMatrix{T}}
    size(C) == size(S) || throw(ArgumentError("the sizes of the center matrix and the " *
                                              "radii matrix should match, but they are $(size(C)) " *
                                              "and $(size(S)) respectively"))

    return IntervalMatrix(map((x, y) -> interval(x - y, x + y), C, S))
end

for op in (:Adjoint, :Bidiagonal, :Diagonal, :Hermitian,
           :SymTridiagonal, :Symmetric, :Transpose, :Tridiagonal)
    @eval LinearAlgebra.$op(A::IntervalMatrix) = IntervalMatrix($op(A.mat))
end

if VERSION >= v"1.3"
    LinearAlgebra.UpperHessenberg(A::IntervalMatrix) = IntervalMatrix(UpperHessenberg(A.mat))
end

@static if vIA >= v"0.22"
    function Base.:(==)(A::IntervalMatrix, B::IntervalMatrix)
        return size(A) == size(B) && all(map((a, b) -> isequal_interval(a, b), A, B))
    end
end

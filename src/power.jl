import Base: copy, get

"""
    IntervalMatrixPower{T}

A wrapper for the matrix power that can be incremented.

### Fields

- `M`  -- the original matrix
- `Mᵏ` -- the current matrix power, i.e., ``M^k``
- `k`  -- the current power index

### Notes

The wrapper should only be accessed using the interface functions.
The internal representation (such as the fields) are subject to future changes.

### Examples

```jldoctest
julia> A = IntervalMatrix([Interval(2.0) Interval(2.0, 3.0); Interval(0.0) Interval(-1.0, 1.0)])
2×2 IntervalMatrix{Float64,Interval{Float64},Array{Interval{Float64},2}}:
 [2, 2]   [2, 3]
 [0, 0]  [-1, 1]

julia> pow = IntervalMatrixPower(A);

julia> increment!(pow)
2×2 IntervalMatrix{Float64,Interval{Float64},Array{Interval{Float64},2}}:
 [4, 4]  [2, 9]
 [0, 0]  [0, 1]

julia> increment(pow)
2×2 IntervalMatrix{Float64,Interval{Float64},Array{Interval{Float64},2}}:
 [8, 8]  [-1, 21]
 [0, 0]   [-1, 1]

julia> get(pow)
2×2 IntervalMatrix{Float64,Interval{Float64},Array{Interval{Float64},2}}:
 [4, 4]  [2, 9]
 [0, 0]  [0, 1]

julia> index(pow)
2

julia> base(pow)
2×2 IntervalMatrix{Float64,Interval{Float64},Array{Interval{Float64},2}}:
 [2, 2]   [2, 3]
 [0, 0]  [-1, 1]

```
"""
mutable struct IntervalMatrixPower{T}
    M::IntervalMatrix{T}
    Mᵏ::IntervalMatrix{T}
    k::Int

    function IntervalMatrixPower(M::IntervalMatrix{T}, Mᵏ::IntervalMatrix{T},
                                 k::Int) where {T}
        @assert k >= 1 "matrix powers must be positive"
        return new{T}(M, Mᵏ, k)
    end
end

function IntervalMatrixPower(M::IntervalMatrix{T}) where {T}
    return IntervalMatrixPower(M, M, 1)
end

function IntervalMatrixPower(M::IntervalMatrix{T}, k::Int) where {T}
    @assert k >= 1 "matrix powers must be positive"
    pow = IntervalMatrixPower(M, M, 1)
    @inbounds for i in 1:(k-1)
        increment!(pow)
    end
    return pow
end

function copy(pow::IntervalMatrixPower)
    return IntervalMatrixPower(pow.M, pow.Mᵏ, pow.k)
end

"""
    increment!(pow::IntervalMatrixPower)

Increment a matrix power in-place (i.e., storing the result in `pow`).

### Input

- `pow` -- wrapper of a matrix power (modified in this function)

### Output

The next matrix power now, reflected in the modified wrapper.
"""
function increment!(pow::IntervalMatrixPower)
    pow.k += 1
    if _isapoweroftwo(pow.k)
        pow.Mᵏ = _eval_poweroftwo(pow.M, pow.k)
    else
        pow.Mᵏ *= pow.M
    end
    return pow.Mᵏ
end

"""
    increment(pow::IntervalMatrixPower)

Increment a matrix power without modifying `pow`.

### Input

- `pow` -- wrapper of a matrix power

### Output

The next matrix power now.
"""
function increment(pow::IntervalMatrixPower)
    return increment!(copy(pow))
end

# checks whether a number is a power of 2
function _isapoweroftwo(k::Integer)
    # see https://stackoverflow.com/a/600306
    return (k & (k - 1)) == 0
end

# compute the (exact) matrix power for a power of two
function _eval_poweroftwo(M::IntervalMatrix, k::Integer)
    Mᵏ = M
    @inbounds for i in 1:Int(log(2, k))
        Mᵏ = square(Mᵏ)
    end
    return Mᵏ
end

"""
    base(pow::IntervalMatrixPower)

Return the original matrix represented by a wrapper of a matrix power.

### Input

- `pow` -- wrapper of a matrix power

### Output

The matrix ``M`` being the basis of the matrix power ``M^k`` represented by the
wrapper.
"""
function base(pow::IntervalMatrixPower)
    return pow.M
end

"""
    get(pow::IntervalMatrixPower)

Return the matrix represented by a wrapper of a matrix power.

### Input

- `pow` -- wrapper of a matrix power

### Output

The matrix power represented by the wrapper.
"""
function get(pow::IntervalMatrixPower)
    return pow.Mᵏ
end

"""
    index(pow::IntervalMatrixPower)

Return the current index of the wrapper of a matrix power.

### Input

- `pow` -- wrapper of a matrix power

### Output

The index `k` of the wrapper representing ``M^k``.
"""
function index(pow::IntervalMatrixPower)
    return pow.k
end

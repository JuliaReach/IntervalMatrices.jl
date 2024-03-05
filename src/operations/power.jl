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
julia> A = IntervalMatrix([interval(2, 2) interval(2, 3); interval(0, 0) interval(-1, 1)])
2×2 IntervalMatrix{Float64, Interval{Float64}, Matrix{Interval{Float64}}}:
 [2.0, 2.0]   [2.0, 3.0]
 [0.0, 0.0]  [-1.0, 1.0]

julia> pow = IntervalMatrixPower(A);

julia> increment!(pow)
2×2 IntervalMatrix{Float64, Interval{Float64}, Matrix{Interval{Float64}}}:
 [4.0, 4.0]  [2.0, 9.0]
 [0.0, 0.0]  [0.0, 1.0]

julia> increment(pow)
2×2 IntervalMatrix{Float64, Interval{Float64}, Matrix{Interval{Float64}}}:
 [8.0, 8.0]  [-1.0, 21.0]
 [0.0, 0.0]  [-1.0, 1.0]

julia> matrix(pow)
2×2 IntervalMatrix{Float64, Interval{Float64}, Matrix{Interval{Float64}}}:
 [4.0, 4.0]  [2.0, 9.0]
 [0.0, 0.0]  [0.0, 1.0]

julia> index(pow)
2

julia> base(pow)
2×2 IntervalMatrix{Float64, Interval{Float64}, Matrix{Interval{Float64}}}:
 [2.0, 2.0]   [2.0, 3.0]
 [0.0, 0.0]  [-1.0, 1.0]

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

function IntervalMatrixPower(M::IntervalMatrix{T}, k::Int;
                             algorithm::String="power") where {T}
    @assert k >= 1 "matrix powers must be positive"
    pow = IntervalMatrixPower(M, M, 1)
    @inbounds for i in 1:(k - 1)
        increment!(pow; algorithm=algorithm)
    end
    return pow
end

function copy(pow::IntervalMatrixPower)
    return IntervalMatrixPower(pow.M, pow.Mᵏ, pow.k)
end

const default_algorithm = "intersect"

"""
    increment!(pow::IntervalMatrixPower; [algorithm=default_algorithm])

Increment a matrix power in-place (i.e., storing the result in `pow`).

### Input

- `pow`       -- wrapper of a matrix power (modified in this function)
- `algorithm` -- (optional; default: `default_algorithm`) algorithm to compute
                 the matrix power; available options:
    * `"multiply"` -- fast computation using `*` from the previous result
    * `"power"` -- recomputation using `^`
    * `"decompose_binary"` -- decompose `k = 2a + b`
    * `"intersect"` -- combination of `"multiply"`/`"power"`/`"decompose_binary"`

### Output

The next matrix power, reflected in the modified wrapper.

### Notes

Independent of `"algorithm"`, if the index is a power of two, we compute the
exact result using squaring.
"""
function increment!(pow::IntervalMatrixPower;
                    algorithm::String=default_algorithm)
    pow.k += 1
    if _isapoweroftwo(pow.k)
        pow.Mᵏ = _eval_poweroftwo(pow)
    elseif algorithm == "multiply"
        pow.Mᵏ = _eval_multiply(pow)
    elseif algorithm == "power"
        pow.Mᵏ = _eval_power(pow)
    elseif algorithm == "decompose_binary"
        pow.Mᵏ = _eval_decompose_binary(pow)
    elseif algorithm == "intersect"
        pow.Mᵏ = _eval_intersect(pow)
    else
        pow.k -= 1
        throw(ArgumentError("algorithm $algorithm is not available; choose " *
                            "from 'multiply', 'power', 'decompose_binary', 'intersect'"))
    end
    return pow.Mᵏ
end

"""
    increment(pow::IntervalMatrixPower; [algorithm=default_algorithm])

Increment a matrix power without modifying `pow`.

### Input

- `pow` -- wrapper of a matrix power
- `algorithm` -- (optional; default: `default_algorithm`) algorithm to compute
                 the matrix power; see [`increment!`](@ref) for available options

### Output

The next matrix power.
"""
function increment(pow::IntervalMatrixPower;
                   algorithm::String=default_algorithm)
    return increment!(copy(pow); algorithm=algorithm)
end

# checks whether a number is a power of 2
function _isapoweroftwo(k::Integer)
    # see https://stackoverflow.com/a/600306
    return (k & (k - 1)) == 0
end

# compute the (exact) matrix power for a power of two
function _eval_poweroftwo(pow::IntervalMatrixPower)
    Mᵏ = pow.M
    @inbounds for i in 1:Int(log(2, pow.k))
        Mᵏ = square(Mᵏ)
    end
    return Mᵏ
end

function _eval_multiply(pow::IntervalMatrixPower)
    return pow.Mᵏ * pow.M
end

function _eval_power(pow::IntervalMatrixPower)
    return pow.M^pow.k
end

function _eval_decompose_binary(pow::IntervalMatrixPower)
    return _eval_decompose_binary_helper(pow.M, pow.k)
end

# decompose k = 2a + b with a, b being positive integers and a being maximal
function _eval_decompose_binary_helper(M, k)
    if k == 1
        Mᵏ = M
    else
        Mᵏ = square(_eval_decompose_binary_helper(M, div(k, 2)))
        if k % 2 == 1
            Mᵏ *= M
        end
    end
    return Mᵏ
end

function _eval_intersect(pow::IntervalMatrixPower)
    return intersect(intersect(_eval_multiply(pow), _eval_power(pow)),
                     _eval_decompose_binary(pow))
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
    matrix(pow::IntervalMatrixPower)

Return the matrix represented by a wrapper of a matrix power.

### Input

- `pow` -- wrapper of a matrix power

### Output

The matrix power represented by the wrapper.
"""
function matrix(pow::IntervalMatrixPower)
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

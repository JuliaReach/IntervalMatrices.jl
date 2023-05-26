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
            @inbounds B[i, j] = _rand_interval(; N=N, rng=rng)
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

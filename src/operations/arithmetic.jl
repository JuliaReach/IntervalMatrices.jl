import Base: +, -, *, /

# =========================
# Addition operations
# =========================

+(M1::IntervalMatrix, M2::IntervalMatrix) = IntervalMatrix(M1.mat + M2.mat)
-(M1::IntervalMatrix, M2::IntervalMatrix) = IntervalMatrix(M1.mat - M2.mat)

+(x::Interval, M::IntervalMatrix) = IntervalMatrix(x .+ M.mat)
+(M::IntervalMatrix, x::Interval) = IntervalMatrix(x .+ M.mat)

+(x::Number, M::IntervalMatrix) = Interval(x) + M
+(M::IntervalMatrix, x::Number) = Interval(x) + M

+(M1::IntervalMatrix, M2::AbstractMatrix) = IntervalMatrix(M1.mat + M2)
+(M1::AbstractMatrix, M2::IntervalMatrix) = IntervalMatrix(M1 + M2.mat)
-(M1::IntervalMatrix, M2::AbstractMatrix) = IntervalMatrix(M1.mat - M2)
-(M1::AbstractMatrix, M2::IntervalMatrix) = IntervalMatrix(M1 - M2.mat)

# =========================
# Multiplication operations
# =========================


*(x::Interval, M::IntervalMatrix) = IntervalMatrix(x .* M.mat)
*(M::IntervalMatrix, x::Interval) = IntervalMatrix(x .* M.mat)

*(x::Number, M::IntervalMatrix) = Interval(x) * M
*(M::IntervalMatrix, x::Number) = Interval(x) * M

/(M::IntervalMatrix, x::Number) = IntervalMatrix(M ./ x)


# =========================
# Multiplication operations
# =========================

"""
    square(A::IntervalMatrix)

Compute the square of an interval matrix.

### Input

- `A` -- interval matrix

### Output

An interval matrix equivalent to `A * A`.

### Algorithm

We follow [1, Section 6].

[1] Kosheleva, Kreinovich, Mayer, Nguyen. Computing the cube of an interval
matrix is NP-hard. SAC 2005.
"""
function square(A::IntervalMatrix)
    B = similar(A)
    n = checksquare(A)

    # case i == j
    @inbounds for j in 1:n
        B[j, j] = pow(A[j, j], 2)
        for k in 1:n
            k == j && continue
            B[j, j] += A[j, k] * A[k, j]
        end
    end

    # case i ≠ j
    @inbounds for j in 1:n
        for i in 1:n
            i == j && continue
            B[i, j] = A[i, j] * (A[j, j] + A[i, i])
            for k in 1:n
                (k == i || k == j) && continue
                B[i, j] += A[i, k] * A[k, j]
            end
        end
    end
    return B
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

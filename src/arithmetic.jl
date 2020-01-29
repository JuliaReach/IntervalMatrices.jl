import Base: +, -, *

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

*(M1::IntervalMatrix, M2::IntervalMatrix) = IntervalMatrix(M1.mat * M2.mat)

*(x::Interval, M::IntervalMatrix) = IntervalMatrix(x .* M.mat)
*(M::IntervalMatrix, x::Interval) = IntervalMatrix(x .* M.mat)

*(x::Number, M::IntervalMatrix) = Interval(x) * M
*(M::IntervalMatrix, x::Number) = Interval(x) * M

*(M1::IntervalMatrix, M2::AbstractMatrix) = IntervalMatrix(M1.mat * M2)
*(M1::AbstractMatrix, M2::IntervalMatrix) = IntervalMatrix(M1 * M2.mat)

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
    B = similar(A.mat)
    n = checksquare(A)
    for j in 1:n
        for i in 1:n
            if i == j
                res = A[i, i]^2
                for k in 1:n
                    if k != i
                        res += A[i, k] * A[k, i]
                    end
                end
                B[i, i] = res
            else
                res = A[i, j] * (A[i, i] + A[j, j])
                for k in 1:n
                    if k != i && k != j
                        res += A[i, k] * A[k, j]
                    end
                end
                B[i, j] = res
            end
        end
    end
    return IntervalMatrix(B)
end

import Base: +, *

# =========================
# Addition operations
# =========================

+(M1::IntervalMatrix, M2::IntervalMatrix) = IntervalMatrix(M1.mat + M2.mat)

+(x::Interval, M::IntervalMatrix) = IntervalMatrix(x .+ M.mat)
+(M::IntervalMatrix, x::Interval) = IntervalMatrix(x .+ M.mat)

+(x::Number, M::IntervalMatrix) = Interval(x) + M
+(M::IntervalMatrix, x::Number) = +(Interval(x), M)

# =========================
# Multiplication operations
# =========================

*(M1::IntervalMatrix, M2::IntervalMatrix) = IntervalMatrix(M1.mat * M2.mat)

function *(x::Interval, A::IntervalMatrix)
    m, n = size(A)
    B = similar(A)
    for j in 1:n
        for i in 1:m
            @inbounds B[i, j] = x * A[i, j]
        end
    end
    return B
end

*(A::IntervalMatrix, x::AbstractInterval) = x * A

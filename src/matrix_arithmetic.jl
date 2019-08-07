import Base: +, *

# =========================
# Addition operations
# =========================

+(M1::IntervalMatrix, M2::IntervalMatrix) = IntervalMatrix(M1.mat + M2.mat)

# =========================
# Multiplication operations
# =========================

*(M1::IntervalMatrix, M2::IntervalMatrix) = IntervalMatrix(M1.mat * M2.mat)

function *(x::ClosedInterval, A::IntervalMatrix)
    m, n = size(A)
    B = Matrix{typeof(x)}(undef, m, n)
    @inbounds for j in 1:n
        for i in 1:m
            B[i, j] = x * A[i, j]
        end
    end
    return IntervalMatrix(B)
end

*(A::IntervalMatrix, x::ClosedInterval) = x * A

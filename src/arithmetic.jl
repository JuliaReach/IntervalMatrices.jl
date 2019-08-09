import Base: +, *

# =========================
# Addition operations
# =========================

+(M1::IntervalMatrix, M2::IntervalMatrix) = IntervalMatrix(M1.mat + M2.mat)

+(x::Interval, M::IntervalMatrix) = IntervalMatrix(x .+ M.mat)
+(M::IntervalMatrix, x::Interval) = IntervalMatrix(x .+ M.mat)

+(x::Number, M::IntervalMatrix) = Interval(x) + M
+(M::IntervalMatrix, x::Number) = Interval(x) + M

# =========================
# Multiplication operations
# =========================

*(M1::IntervalMatrix, M2::IntervalMatrix) = IntervalMatrix(M1.mat * M2.mat)

*(x::Interval, M::IntervalMatrix) = IntervalMatrix(x .* M.mat)
*(M::IntervalMatrix, x::Interval) = IntervalMatrix(x .* M.mat)

# =========================
# Multiplication operations
# =========================
import Base.*

function *(λ::Number, x::ClosedInterval)
    return ClosedInterval(λ*x.left, λ*x.right)
end

function *(x::ClosedInterval, λ::Number)
    return ClosedInterval(λ*x.left, λ*x.right)
end

function *(x::ClosedInterval, y::ClosedInterval)
    z = (x.left * y.left, x.left * y.right, x.right * y.left, x.right * y.right)
    return ClosedInterval(minimum(z), maximum(z))
end

# ===================
# Addition operations
# ===================
import Base.+

function +(λ::Number, x::ClosedInterval)
    return ClosedInterval(λ + x.left, λ + x.right)
end

function +(x::ClosedInterval, λ::Number)
    return ClosedInterval(λ + x.left, λ + x.right)
end

function +(x::ClosedInterval, y::ClosedInterval)
    return ClosedInterval(x.left + y.left, x.right + y.right)
end

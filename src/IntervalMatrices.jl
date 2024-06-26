module IntervalMatrices

using LinearAlgebra
using LinearAlgebra: checksquare

using Random: AbstractRNG, GLOBAL_RNG, seed!
import Random: rand

import Base: copy, get,
             +, -, *, /, \,
             size, IndexStyle, getindex, setindex!,
             similar, ∈, ⊆, ∩, ∪, real, imag

using Reexport
@reexport using IntervalArithmetic
import IntervalArithmetic: inf, sup, mid, diam, radius, hull
@static if VERSION >= v"1.9"
    vIA = pkgversion(IntervalArithmetic)
else
    import PkgVersion
    vIA = PkgVersion.Version(IntervalArithmetic)
end
if vIA >= v"0.21"
    # IntervalArithmetic v0.21 removed `convert`
    Base.convert(::Type{Interval{T}}, x::Number) where {T} = interval(T(x))
    Base.convert(::Type{Interval{T}}, x::Interval{T}) where {T} = x
    function Base.convert(::Type{Interval{T1}}, x::Interval{T2}) where {T1,T2}
        return interval(T1(inf(x)), T1(sup(x)))
    end
else
    # COV_EXCL_START
    # IntervalArithmetic v0.21 requires `interval`, but prior versions did not
    # offer this method for `Complex` inputs
    IntervalArithmetic.interval(a::Complex) = complex(interval(real(a)), interval(imag(a)))
    # COV_EXCL_STOP
end
if vIA >= v"0.22"
    import Base: intersect
    export ±, midpoint_radius

    function ±(x::Number, y::Number)
        return x + interval(-y, y)
    end

    function Base.:(==)(A::AbstractMatrix{<:Interval}, B::AbstractMatrix{<:Interval})
        return size(A) == size(B) && all(map((a, b) -> isequal_interval(a, b), A, B))
    end

    function intersect(a::Interval{T}, b::Interval{T}) where {T}
        lo = max(inf(a), inf(b))
        hi = min(sup(a), sup(b))
        if lo > hi
            return emptyinterval(T)
        end
        return interval(lo, hi)
    end
    intersect(a::Interval{T}, b::Interval{S}) where {T,S} = intersect(promote(a, b)...)

    Base.in(x::Number, y::Interval) = inf(y) <= x <= sup(y)
else
    import IntervalArithmetic: ±, midpoint_radius

    issubset_interval(x::Interval, y::Interval) = issubset(x, y)
end

# ================================
# Type defining an interval matrix
# ================================
include("matrix.jl")
include("affine.jl")

# =================================
# Operations for interval matrices
# =================================
include("operations/arithmetic.jl")
include("operations/mult.jl")
include("operations/norm.jl")
include("operations/numops.jl")
include("operations/random.jl")
include("operations/setops.jl")

# ==================================================
# Methods to compute the power of an interval matrix
# ==================================================
include("operations/power.jl")

# =======================================================
# Methods to handle the exponential of an interval matrix
# =======================================================
include("exponential.jl")

# ==============================================================
# Methods to compute a correction matrix from an interval matrix
# ==============================================================
include("correction_matrices.jl")

# ========
# Exports
# ========
export AbstractIntervalMatrix,
       IntervalMatrix,
       set_multiplication_mode, get_multiplication_mode

export inf, sup,
       rand, sample,
       opnorm,
       diam_norm,
       scale, scale!,
       square,
       quadratic_expansion,
       exp_overapproximation,
       exp_underapproximation

export IntervalMatrixPower,
       increment!, increment,
       base,
       index

export AffineIntervalMatrix,
       AffineIntervalMatrix1

set_multiplication_mode(config[:multiplication])

end # module

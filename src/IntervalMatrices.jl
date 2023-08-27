__precompile__(true)
module IntervalMatrices

using LinearAlgebra
using LinearAlgebra: checksquare

import Random
using Random: AbstractRNG, GLOBAL_RNG, seed!

using Reexport
@reexport using IntervalArithmetic
@static if VERSION >= v"1.9"
    vIA = pkgversion(IntervalArithmetic)
else
    using PkgVersion
    vIA = PkgVersion.Version(IntervalArithmetic)
end
if vIA >= v"0.21"
    # IntervalArithmetic v0.21 removed convert
    Base.convert(::Type{Interval{T}}, x::Number) where {T} = interval(T(x))
    Base.convert(::Type{Interval{T}}, x::Interval{T}) where {T} = x
    Base.convert(::Type{Interval{T}}, x::Interval) where {T} = interval(T(inf(x)), T(sup(x)))
else
    # IntervalArithmetic v0.21 requires interval, but prior versions did not offer this method
    IntervalArithmetic.interval(a::Complex) = complex(interval(real(a)), interval(imag(a)))
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

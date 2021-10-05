__precompile__(true)
module IntervalMatrices

using LinearAlgebra

import Random
using Random: AbstractRNG, GLOBAL_RNG, seed!

using Reexport
@reexport using IntervalArithmetic

# ================================
# Type defining an interval matrix
# ================================
include("matrix.jl")

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


set_multiplication_mode(config[:multiplication])


end # module

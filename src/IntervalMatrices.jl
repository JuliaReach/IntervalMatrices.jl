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
# Arithmetic for interval matrices
# =================================
include("arithmetic.jl")

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
       IntervalMatrix

export inf, sup,
       rand, sample,
       opnorm,
       quadratic_expansion,
       expm_overapproximation,
       expm_underapproximation

end # module

__precompile__(true)

"""
Main module for `IntervalMatrices.jl` -- a Julia package to work with
matrices that have uncertain parameters.
"""
module IntervalMatrices

using LinearAlgebra

import Random
using Random: AbstractRNG, GLOBAL_RNG, seed!

using Reexport
@reexport using IntervalSets

# ================================
# Arithmetic for interval sets
# ================================
include("arithmetic.jl")

# ================================
# Type defining an interval matrix
# ================================
include("matrix.jl")

# =======================================================
# Methods to handle the exponential of an interval matrix
# =======================================================
include("exponential.jl")

# ========
# Exports
# ========
export AbstractIntervalMatrix,
       IntervalMatrix

export opnorm,
       left,
       right,
       quadratic_expansion,
       expm_overapproximation,
       expm_underapproximation

end # module

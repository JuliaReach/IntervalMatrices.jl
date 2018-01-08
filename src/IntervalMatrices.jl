__precompile__(true)

"""
Main module for `IntervalMatrices.jl` -- a Julia package to work with
matrices that have uncertain parameters.
"""
module IntervalMatrices

using Reexport
@reexport using IntervalSets

include("arithmetic.jl")
include("matrix.jl")
include("exponential.jl")

# ======
# Types
# ======

export AbstractIntervalMatrix,
       IntervalMatrix,
       IntervalMatrixExponential

# ========
# Methods
# ========
export left,
       right,
       quadratic_expansion,
       expm_overapproximation,
       expm_underapproximation

end # module

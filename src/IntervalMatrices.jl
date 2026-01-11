module IntervalMatrices

import Base: copy, rand,
             +, -, *, /, \,
             size, IndexStyle, getindex, setindex!,
             similar, ∈, ⊆, ∩, ∪, real, imag

import LinearAlgebra
using LinearAlgebra: Diagonal, I, UniformScaling, checksquare, opnorm
using Random: AbstractRNG, GLOBAL_RNG
using Reexport: @reexport

# =================================
# Interface with IntervalArithmetic
# =================================
include("init_IntervalArithmetic.jl")

# ================================
# Type defining an interval matrix
# ================================
include("matrix.jl")
include("affine.jl")

# ================================
# Operations for interval matrices
# ================================
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

# =======
# Exports
# =======
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
       matrix,
       base,
       index

export AffineIntervalMatrix,
       AffineIntervalMatrix1

set_multiplication_mode(config[:multiplication])

function __init__()
    @static if vIA >= v"0.22"
        # remove interval decorations
        setdisplay(; decorations=false, ng_flag=false)
    end
end

end  # module

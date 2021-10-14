using IntervalMatrices, Test, LinearAlgebra, SparseArrays

using IntervalMatrices: _truncated_exponential_series,
                        horner,
                        scale_and_square,
                        correction_hull,
                        input_correction
include("models.jl")

include("constructors.jl")
include("arithmetic.jl")
include("setops.jl")
include("exponential.jl")
include("affine.jl")

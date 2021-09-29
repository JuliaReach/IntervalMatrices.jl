using IntervalMatrices, Test, LinearAlgebra

using IntervalMatrices: _truncated_exponential_series,
                        horner,
                        scale_and_square,
                        correction_hull,
                        input_correction

include("constructors.jl")
include("arithmetic.jl")
include("setops.jl")
include("exponential.jl")

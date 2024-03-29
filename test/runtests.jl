using IntervalMatrices, Test, LinearAlgebra, SparseArrays
import PkgVersion

using IntervalMatrices: _truncated_exponential_series,
                        horner,
                        scale_and_square,
                        correction_hull,
                        input_correction

@static if PkgVersion.Version(IntervalMatrices.IntervalArithmetic) >= v"0.22"
    # equality test for convenience
    Base.:(==)(x::Interval, y::Interval) = isequal_interval(x, y)
end

include("models.jl")

include("constructors.jl")
include("arithmetic.jl")
include("setops.jl")
include("exponential.jl")
include("affine.jl")

include("Aqua.jl")

@reexport using IntervalArithmetic
import IntervalArithmetic: inf, sup, mid, diam, radius, hull

@static if VERSION >= v"1.9"
    vIA = pkgversion(IntervalArithmetic)
else
    import PkgVersion
    vIA = PkgVersion.Version(IntervalArithmetic)
end

@static if vIA >= v"0.21.2"
    import IntervalArithmetic: midradius
else
    export midradius
end

@static if vIA >= v"0.22"
    import Base: intersect
    export ±  # not exported by IntervalArithmetic anymore
else  # vIA < v"0.22"
    # COV_EXCL_START
    import IntervalArithmetic: ±

    issubset_interval(x::Interval, y::Interval) = issubset(x, y)

    in_interval(x::Number, y::Interval) = inf(y) <= x <= sup(y)

    intersect_interval(a::Interval, b::Interval) = intersect(a, b)

    if vIA >= v"0.21"
        # `convert` was temporarily removed in IntervalArithmetic v0.21 until v0.22
        Base.convert(::Type{Interval{T}}, x::Number) where {T} = interval(T(x))
        Base.convert(::Type{Interval{T}}, x::Interval{T}) where {T} = x
        function Base.convert(::Type{Interval{T1}}, x::Interval{T2}) where {T1,T2}
            return interval(T1(inf(x)), T1(sup(x)))
        end
    else  # vIA < v"0.21"
        # IntervalArithmetic v0.21 requires `interval`, but prior versions did not
        # offer this method for `Complex` inputs
        IntervalArithmetic.interval(a::Complex) = complex(interval(real(a)), interval(imag(a)))
    end
    # COV_EXCL_STOP
end

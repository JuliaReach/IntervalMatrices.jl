# IntervalMatrices.jl

`IntervalMatrices` is a [Julia](http://julialang.org) package to work with
matrices that have interval coefficients.

## Features

Here is a quick summary of the available functionality.
See the section [Library Outline](@ref) below for details.

- An [`IntervalMatrix`](@ref) type that wraps a two-dimensional array whose components
  are intervals.
- Arithmetic between scalars, intervals and interval matrices.
- Quadratic expansions using single-use expressions (SUE).
- Interval matrix exponential: underapproximation and overapproximation routines.
- Utility functions such as: operator norm, random sampling, splitting and
  containment check.

An application of interval matrices is to find the set of states reachable by
a dynamical system whose coefficients are uncertain. The library
[ReachabilityAnalysis.jl](http://github.com/JuliaReach/ReachabilityAnalysis.jl)
implements algorithms that use interval matrices [AlthoffKS11, Liou66](@cite).

## Installing

Depending on your needs, choose an appropriate command from the following list
and enter it in Julia's REPL.

To [install the latest release version](https://julialang.github.io/Pkg.jl/v1/managing-packages/#Adding-registered-packages-1):

```
pkg> add IntervalMatrices
```

To install the latest development version:

```
pkg> add IntervalMatrices#master
```

To [clone the package for development](https://julialang.github.io/Pkg.jl/v1/managing-packages/#Developing-packages-1):

```
pkg> dev IntervalMatrices
```

## Quickstart

An *interval matrix* is a matrix whose coefficients are intervals. For instance,

```jldoctest quickstart
julia> using IntervalMatrices

julia> A = IntervalMatrix([interval(0, 1) interval(1, 2); interval(2, 3) interval(-4, -2)])
2×2 IntervalMatrix{Float64, Interval{Float64}, Matrix{Interval{Float64}}}:
 [0.0, 1.0]   [1.0, 2.0]
 [2.0, 3.0]  [-4.0, -2.0]
```
defines an interval matrix ``A``. The type of each coefficient in ``A`` is an interval,
e.g. its coefficient in position ``(1, 1)`` is the interval ``[0, 1]`` over double-precision
floating-point numbers:

```jldoctest quickstart
julia> A[1, 1]
[0.0, 1.0]

julia> typeof(A[1, 1])
Interval{Float64}
```
This library uses the interval arithmetic package
[IntervalArithmetic.jl](https://github.com/JuliaIntervals/IntervalArithmetic.jl)
to deal with interval computations. For instance, one can compute a multiple
of ``A``,

```jldoctest quickstart
julia> 2A
2×2 IntervalMatrix{Float64, Interval{Float64}, Matrix{Interval{Float64}}}:
 [0.0, 2.0]   [2.0, 4.0]
 [4.0, 6.0]  [-8.0, -4.0]
```
Or an interval multiple of ``A``,

```jldoctest quickstart
julia> interval(-1, 1) * A
2×2 IntervalMatrix{Float64, Interval{Float64}, Matrix{Interval{Float64}}}:
 [-1.0, 1.0]  [-2.0, 2.0]
 [-3.0, 3.0]  [-4.0, 4.0]
```

Or compute the square of ``A``,
```jldoctest quickstart
julia> square(A)
2×2 IntervalMatrix{Float64, Interval{Float64}, Matrix{Interval{Float64}}}:
   [2.0, 7.0]   [-8.0, -1.0]
 [-12.0, -2.0]   [6.0, 22.0]
```
In these cases, the rules of interval arithmetic are used; see the Wikipedia page
on [interval arithmetic](https://en.wikipedia.org/wiki/Interval_arithmetic) for the
relevant definitions and algebraic rules that apply.

However, the straightforward application of the rules of interval arithmetic does
not always give the exact result: in general it only gives an overapproximation
[AlthoffSB08, KoshelevaKMN05](@cite). To illustrate, suppose that we
are interested in the quadratic term ``At + \frac{1}{2}A^2 t^2``, which
corresponds to the Taylor-series expansion at order two of ``e^{At} - I``.
Then, at ``t = 1.0``,

```jldoctest quickstart
julia> A + 1/2 * A^2
2×2 IntervalMatrix{Float64, Interval{Float64}, Matrix{Interval{Float64}}}:
  [0.5, 4.5]   [-3.25, 2.5]
 [-4.25, 3.0]  [-2.25, 9.0]
```
However, that result is not tight. The computation can be performed exactly via
single-use expressions implemented in this library:

```jldoctest quickstart
julia> quadratic_expansion(A, 1.0, 0.5)
2×2 IntervalMatrix{Float64, Interval{Float64}, Matrix{Interval{Float64}}}:
  [1.0, 4.5]  [-2.0, 1.0]
 [-3.0, 1.5]   [1.0, 7.0]
```
We now obtain an interval matrix that is strictly included in the one obtained from
the naive multiplication.

An overapproximation and an underapproximation method at a given order for
``e^{At}``, where ``A`` is an interval matrix, are also available. See the [Methods](@ref)
section for details.

## Library Outline

Explore the types and methods defined in this library by following the links below,
or use the search bar in the left to look for a specific keyword in the documentation.

```@contents
Pages = [
    "lib/types.md",
    "lib/methods.md"
]
Depth = 2
```

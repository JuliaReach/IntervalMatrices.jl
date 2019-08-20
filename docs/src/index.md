# IntervalMatrices.jl

`IntervalMatrices` is a [Julia](http://julialang.org) package to work with
matrices that have uncertain parameters.

## Features

Here is a quick summary of the available functionality.
See the section `Library Outline` below for details.

- An `IntervalMatrix` type that wraps a two-dimensional array whose components
  are intervals.
- Arithmetic between scalars, intervals and interval matrices.
- Quadratic expansions using single-use expressions (SUE).
- Interval matrix exponential: underapproximation and overapproximation routines.
- Utility functions such as: operator norm, random sampling, splitting and
  containment check.

## Quickstart

An *interval matrix* is a matrix whose coefficients are intervals. For instance,

```jldoctest quickstart
julia> using IntervalMatrices

julia> A = IntervalMatrix([0..1 1..2; 2..3 -4.. -2])
2×2 IntervalMatrix{Float64,Interval{Float64},Array{Interval{Float64},2}}:
 [0, 1]    [1, 2]
 [2, 3]  [-4, -2]
```
defines an interval matrix $A$. The type of each coefficient in $A$ is an interval,
e.g. its coefficient in position $(1, 1)$ is the interval $[0, 1]$ over double-precision
floating-point numbers:

```jldoctest quickstart
julia> A[1, 1]
[0, 1]

julia> typeof(A[1, 1])
Interval{Float64}
```
This library uses the interval arithmetic package
[IntervalArithmetic.jl](https://github.com/JuliaIntervals/IntervalArithmetic.jl)
to deal with interval computations. For instance, one can compute a multiple
of $A$

```jldoctest quickstart
julia> 2A
2×2 Array{Interval{Float64},2}:
 [0, 2]    [2, 4]
 [4, 6]  [-8, -4]
```
Or an interval multiple of $A$,

```jldoctest quickstart
julia> (-1.0..1.0) * A
2×2 IntervalMatrix{Float64,Interval{Float64},Array{Interval{Float64},2}}:
 [-1, 1]  [-2, 2]
 [-3, 3]  [-4, 4]
```

Or compute the square of $A$,
```jldoctest quickstart
julia> A*A
2×2 IntervalMatrix{Float64,Interval{Float64},Array{Interval{Float64},2}}:
    [2, 7]  [-8, 0]
 [-12, -1]  [6, 22]
```
In these cases, the rules of interval arithmetic are used; see the wikipedia page
on [interval arithmetic](https://en.wikipedia.org/wiki/Interval_arithmetic) for the
relevant definitions and algebraic rules that apply.

However, the straightforward application of the rules of interval arithmetic does
not always give the exact result; in general it only gives an overapproximation [1, 2].
To illustrate, suppose that we are interested in the quadratic term
$At + \frac{1}{2}A^2 t^2$, which corresponds to the Taylor-series expansion at order two of
$e^{At} - I$. Then, at $t = 1.0$,

```jldoctest quickstart
julia> A + 1/2 * A^2
2×2 Array{Interval{Float64},2}:
  [1, 4.5]  [-3, 2]
 [-4, 2.5]  [-1, 9]
```
This computation can be performed exactly via single-use expressions implemented
in `IntervalMatrices.jl`, obtaining an interval matrix that is strictly included
in the previous result:

```jldoctest quickstart
julia> quadratic_expansion(A, 1.0)
2×2 Array{Interval{Float64},2}:
  [1, 4.5]  [-2, 1]
 [-3, 1.5]   [1, 7]
```
An overapproximation and an underapproximation method at a given order for
$e^{At}$, where $A$ is an interval matrix, are also available. See the `Methods`
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

## References

- [1] Althoff, Matthias, Olaf Stursberg, and Martin Buss. *Reachability analysis of nonlinear systems
      with uncertain parameters using conservative linearization.*
      2008 47th IEEE Conference on Decision and Control. IEEE, 2008.

- [2] Kosheleva, Olga, et al. *Computing the cube of an interval matrix is NP-hard.*
      Proceedings of the 2005 ACM symposium on Applied computing. ACM, 2005.

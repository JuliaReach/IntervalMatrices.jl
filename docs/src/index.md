# IntervalMatrices.jl

```@meta
DocTestFilters = [r"[0-9\.]+ seconds \(.*\)"]
```

`IntervalMatrices` is a [Julia](http://julialang.org) package to work with
matrices that have uncertain parameters.

## Features

- An `IntervalMatrix` type.
- Arithmetics between scalars, closed intervals and interval matrices.
- Overapproximation and underapproximation routines for the exponential of an
  interval matrix at different orders.

## Library Outline

```@docs
IntervalMatrices
```

```@contents
Pages = [
    "lib/types.md",
    "lib/methods.md"
]
Depth = 2
```

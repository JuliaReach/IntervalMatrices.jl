# Methods

This section describes systems methods implemented in `IntervalMatrices.jl`.

```@contents
Pages = ["methods.md"]
Depth = 3
```

```@meta
CurrentModule = IntervalMatrices
```

## Common functions

```@docs
inf
sup
mid
diam
rand
sample
split
∈
square
scale
scale!
⊆
∩
```

## Matrix power

```@docs
increment!
increment
get
base
index
```

## Exponentiation

```@docs
exp_overapproximation
horner
scale_and_square
exp_underapproximation
```

## Finite expansions

```@docs
quadratic_expansion
```

## Correction terms

```@docs
correction_hull
input_correction
```

## Norms

```@docs
opnorm
diam_norm
```

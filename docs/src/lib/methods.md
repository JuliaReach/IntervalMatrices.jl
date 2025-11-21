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
radius
midradius
rand
sample
∈
±
⊆
∩
∪
hull
```

## Arithmetic

```@docs
square
scale
scale!
set_multiplication_mode
```

## Matrix power

```@docs
increment!
increment
matrix
base
index
```

## Matrix exponential

### Algorithms

```@docs
IntervalMatrices.Horner
IntervalMatrices.ScaleAndSquare
IntervalMatrices.TaylorOverapproximation
IntervalMatrices.TaylorUnderapproximation
```

### Implementations

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

"""
    opnorm(A::IntervalMatrix, p::Real=Inf)

The matrix norm of an interval matrix.

### Input

- `A` -- interval matrix
- `p` -- (optional, default: `Inf`) the class of `p`-norm

### Notes

The matrix ``p``-norm of an interval matrix ``A`` is defined as

```math
    ‖A‖_p := ‖\\max(|\\text{inf}(A)|, |\\text{sup}(A)|)‖_p
```

where ``\\max`` and ``|·|`` are taken elementwise.
"""
function LinearAlgebra.opnorm(A::IntervalMatrix, p::Real=Inf)
    if p == Inf
        return _opnorm_inf(A)
    elseif p == 1
        return _opnorm_1(A)
    else
        error("the interval matrix norm for this value of p=$p is not implemented")
    end
end

# The operator norm in the infinity norm corresponds to the
# maximum absolute value row-sum.
function _opnorm_inf(A::IntervalMatrix{T}) where {T}
    m, n = size(A)
    res = zero(T)
    @inbounds @simd for i in 1:m
        acc = zero(T)
        for j in 1:n
            x = A[i, j]
            acc += max(abs(inf(x)), abs(sup(x)))
        end
        if acc > res
            res = acc
        end
    end
    return res
end

# The operator norm in the 1-norm corresponds to the
# maximum absolute value column-sum.
function _opnorm_1(A::IntervalMatrix{T}) where {T}
    m, n = size(A)
    res = zero(T)
    @inbounds @simd for j in 1:n
        acc = zero(T)
        for i in 1:m
            x = A[i, j]
            acc += max(abs(inf(x)), abs(sup(x)))
        end
        if acc > res
            res = acc
        end
    end
    return res
end

"""
    diam_norm(A::IntervalMatrix, p=Inf)

Return the diameter norm of the interval matrix.

### Input

- `A` -- interval matrix
- `p` -- (optional, default: `Inf`) the `p`-norm used; valid options are:
         `1`, `2`, `Inf`

### Output

The operator norm, in the `p`-norm, of the scalar matrix obtained
by taking the element-wise `diam` function, where `diam(x) := sup(x) - inf(x)`
for an interval `x`.

### Notes

This function gives a measure of the *width* of the interval matrix.
"""
function diam_norm(A::IntervalMatrix, p=Inf)
    return opnorm(diam.(A), p)
end

"""
    ∈(M::AbstractMatrix, A::AbstractIntervalMatrix)

Check whether a concrete matrix is an instance of an interval matrix.

### Input

- `M` -- concrete matrix
- `A` -- interval matrix

### Output

`true` iff `M` is an instance of `A`

### Algorithm

We check for each entry in `M` whether it belongs to the corresponding interval
in `A`.
"""
function ∈(M::AbstractMatrix, A::AbstractIntervalMatrix)
    @assert size(M) == size(A) "incompatible matrix sizes (M: $(size(M)), A: " *
                               "$(size(A)))"

    m, n = size(A)
    @inbounds for j in 1:n
        for i in 1:m
            if !in_interval(M[i, j], A[i, j])
                return false
            end
        end
    end
    return true
end

"""
    ⊆(A::AbstractIntervalMatrix, B::AbstractIntervalMatrix)

Check whether an interval matrix is contained in another interval matrix.

### Input

- `A` -- interval matrix
- `B` -- interval matrix

### Output

`true` iff `A[i, j] ⊆ B[i, j]` for all `i, j`.
"""
function ⊆(A::AbstractIntervalMatrix, B::AbstractIntervalMatrix)
    @assert size(A) == size(B) "incompatible matrix sizes $(size(A)) and " *
                               "$(size(B))"

    m, n = size(A)
    @inbounds for j in 1:n, i in 1:m
        if !issubset_interval(A[i, j], B[i, j])
            return false
        end
    end
    return true
end

"""
    ∩(A::IntervalMatrix, B::IntervalMatrix)

Intersect two interval matrices.

### Input

- `A` -- interval matrix
- `B` -- interval matrix (of the same shape as `A`)

### Output

A new matrix `C` of the same shape as `A` such that
`C[i, j] = A[i, j] ∩ B[i, j]` for each `i` and `j`.
"""
function ∩(A::IntervalMatrix, B::IntervalMatrix)
    @assert size(A) == size(B) "incompatible matrix sizes (A: $(size(A)), B: " *
                               "$(size(B)))"

    return IntervalMatrix(map((x, y) -> intersect_interval(x, y), A, B))
end

"""
    hull(A::IntervalMatrix, B::IntervalMatrix)

Finds the interval hull of two interval matrices. This is equivalent to [`∪`](@ref).

### Input

- `A` -- interval matrix
- `B` -- interval matrix (of the same shape as `A`)

### Output

A new matrix `C` of the same shape as `A` such that
`C[i, j] = hull(A[i, j], B[i, j])` for each `i` and `j`.
"""
function hull(A::IntervalMatrix, B::IntervalMatrix)
    @assert size(A) == size(B) "incompatible matrix sizes (A: $(size(A)), B: " *
                               "$(size(B)))"

    return IntervalMatrix(map((x, y) -> hull(x, y), A, B))
end

"""
    ∪(A::IntervalMatrix, B::IntervalMatrix)

Finds the interval union (hull) of two interval matrices.
This is equivalent to [`hull`](@ref).

### Input

- `A` -- interval matrix
- `B` -- interval matrix (of the same shape as `A`)

### Output

A new matrix `C` of the same shape as `A` such that
`C[i, j] = A[i, j] ∪ B[i, j]` for each `i` and `j`.
"""
∪(A::IntervalMatrix, B::IntervalMatrix) = hull(A, B)

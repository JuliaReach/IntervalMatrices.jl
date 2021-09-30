"""
    inf(A::IntervalMatrix{T}) where {T}

Return the infimum of an interval matrix `A`, which corresponds to taking the
element-wise infimum of `A`.

### Input

- `A` -- interval matrix

### Output

A scalar matrix whose coefficients are the infima of each element in `A`.
"""
function inf(A::IntervalMatrix{T}) where {T}
    return map(inf, A.mat)
end

"""
    sup(A::IntervalMatrix{T}) where {T}

Return the supremum of an interval matrix `A`, which corresponds to taking the
element-wise supremum of `A`.

### Input

- `A` -- interval matrix

### Output

A scalar matrix whose coefficients are the suprema of each element in `A`.
"""
function sup(A::IntervalMatrix{T}) where {T}
    return map(sup, A.mat)
end

"""
    mid(A::IntervalMatrix{T}) where {T}

Return the midpoint of an interval matrix `A`, which corresponds to taking the
element-wise midpoint of `A`.

### Input

- `A` -- interval matrix

### Output

A scalar matrix whose coefficients are the midpoints of each element in `A`.
"""
function mid(A::IntervalMatrix{T}) where {T}
    return map(mid, A.mat)
end

"""
    radius(A::IntervalMatrix{T}) where {T}

Return the radius of an interval matrix `A`, which corresponds to taking the
element-wise radius of `A`.

### Input

- `A` -- interval matrix

### Output

A scalar matrix whose coefficients are the radii of each element in `A`.
"""
function radius(A::IntervalMatrix{T}) where {T}
    return map(radius, A.mat)
end

"""
    diam(A::IntervalMatrix{T}) where {T}

Return a matrix whose entries describe the diameters of the intervals.

### Input

- `A` -- interval matrix

### Output

A matrix `B` of the same shape as `A` such that `B[i, j] == diam(A[i, j])` for
each `i` and `j`.
"""
function diam(A::IntervalMatrix{T}) where {T}
    return map(diam, A.mat)
end


"""
    midpoint_radius(A::IntervalMatrix{T}) where {T}

Split an interval matrix ``A`` into two scalar matrices ``C`` and ``S``
such that ``A = C + [-S, S]``.

### Input

- `A` -- interval matrix

### Output

A pair `(C, S)` such that the entries of `C` are the central points and the
entries of `S` are the (nonnegative) radii of the intervals in `A`.
"""
midpoint_radius(A::IntervalMatrix) = (mid(A), radius(A))

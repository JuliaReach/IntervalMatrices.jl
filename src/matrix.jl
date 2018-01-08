abstract type AbstractIntervalMatrix{IT} <: AbstractMatrix{IT} end

struct IntervalMatrix{T, IT <: AbstractInterval{T}} <: AbstractIntervalMatrix{IT}
    mat::AbstractMatrix{IT}
end
IntervalMatrix{T, IT <: AbstractInterval{T}}(x::AbstractMatrix{IT}) = IntervalMatrix{T, IT}(x)

import Base:size, IndexStyle, getindex, setindex!, +, *

IndexStyle(::Type{<:IntervalMatrix}) = IndexLinear()
size(M::IntervalMatrix) = size(M.mat)
getindex(M::IntervalMatrix, i::Int) = getindex(M.mat, i)
setindex!(M::IntervalMatrix, X, inds...) = setindex!(M.mat, X, inds...)

+(M1::IntervalMatrix, M2::IntervalMatrix) = IntervalMatrix(M1.mat + M2.mat)
*(M1::IntervalMatrix, M2::IntervalMatrix) = IntervalMatrix(M1.mat * M2.mat)

"""
    norm(A::IntervalMatrix, p::Real=Inf)
"""
function Base.LinAlg.norm(A::IntervalMatrix, p::Real=Inf)
    if p == Inf
        return norm(max.(abs.(left(A)), abs.(right(A))), Inf)
    else
        error("the interval matrix norm for this value of p=$p is not implemented")
    end
end

"""
    left(A::IntervalMatrix)
"""
function left(A::IntervalMatrix)
    n = size(A, 1)
    return hcat([[A[i, j].left for j in 1:n] for i in 1:n]...)'
end

"""
    right(A::IntervalMatrix)
"""
function right(A::IntervalMatrix)
    n = size(A, 1)
    return hcat([[A[i, j].right for j in 1:n] for i in 1:n]...)'
end

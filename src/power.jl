import Base: ^, show

# helper types representing factors

abstract type AbstractFactor{T<:Integer} end

struct OneFactor{T} <: AbstractFactor{T} end

show(io::IO, ::OneFactor) = print(io, 1)

int(::OneFactor{T}) where {T} = one(T)

struct DoubledFactor{T} <: AbstractFactor{T}
    a::T
end

show(io::IO, df::DoubledFactor) = print(io, df.a, "·2")

function int(df::DoubledFactor)
    return df.a * 2
end

# binary decomposition

function binary_decomposition(k::T) where {T<:Integer}
    factor2source = Dict{T, Vector{AbstractFactor{T}}}()
    binary_decomposition!(factor2source, k)
    return factor2source
end

function binary_decomposition!(factor2source, k::T) where {T<:Integer}
    factors = Vector{AbstractFactor{T}}()
    @assert !haskey(factor2source, k)
    factor2source[k] = factors
    if k == 1
        push!(factors, OneFactor{T}())
        return
    else
        b = div(k, 2)
        push!(factors, DoubledFactor(b))
        if !haskey(factor2source, b)
            binary_decomposition!(factor2source, b)
        end
        if mod(k, 2) == 1
            # k = 2 * b + 1
            push!(factors, OneFactor{T}())
        end
    end
    return factor2source
end

# matrix power

function matpow!(factor2matrix, factor2source, k::Integer)
    if haskey(factor2matrix, k)
        return factor2matrix[k]
    end
    A = nothing
    for a in factor2source[k]
        nₐ = int(a)
        if !haskey(factor2matrix, nₐ)
            @assert a isa DoubledFactor
            # B = C² where C needs to be computed recursively
            C = matpow!(factor2matrix, factor2source, a.a)
            B = square(C)
            @assert !haskey(factor2matrix, nₐ)
            factor2matrix[nₐ] = B
        else
            B = factor2matrix[nₐ]
        end
        if A == nothing
            A = B
        else
            A *= B
        end
    end
    factor2matrix[k] = A
    return A
end

function ^(M::IntervalMatrix{T1}, k::T2) where {T1, T2<:Integer}
    factor2source = binary_decomposition(k)
    factor2source[1] = [OneFactor{T2}()]
    factor2matrix = Dict{T2, IntervalMatrix{T1}}(1 => M)
    return matpow!(factor2matrix, factor2source, k)
end

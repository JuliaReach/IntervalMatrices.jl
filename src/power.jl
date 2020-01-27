import Base: ^

function ^(M::IntervalMatrix{T}, k::Integer) where {T}
    n = checksquare(M)

    # create symbolic matrix
    Ms = Matrix{SymEngine.Basic}(undef, n, n)
    dict = Dict{SymEngine.Basic, Interval{T}}()
    for i in 1:n, j in 1:n
        var = SymEngine.symbols("M_$(i)_$(j)")
        Ms[i, j] = var
        dict[var] = M[i, j]
    end

    # compute symbolic matrix power
    Msᵏ = Ms^k

    # substitute intervals in symbolic result
    Mᵏ = similar(M)
    for i in 1:n, j in 1:n
        println(Msᵏ[i, j])
        println(dict)
        substituted = SymEngine.subs(Msᵏ[i, j], dict)
        Mᵏ[i, j] = SymEngine.lambdify(substituted)
    end
    return Mᵏ
end

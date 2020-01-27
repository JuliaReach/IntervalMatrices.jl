import Base: ^

function ^(M::IntervalMatrix{T}, k::Integer) where {T}
    n = checksquare(M)
#     Ms = [SymEngine.symbols("M_$(i)_$(j)") for i in 1:n, j in 1:n]
    Ms = Matrix{SymEngine.Basic}(undef, n, n)
    dict = Dict{SymEngine.Basic, T}()
    for i in 1:n, j in 1:n
        var = SymEngine.symbols("M_$(i)_$(j)")
        Ms[i, j] = var
        dict[var] = M[i, j]
    end
    println(Ms)
    Msᵏ = Ms^k
    Mᵏ = subs(Msᵏ, dict)
    return Mᵏ
end

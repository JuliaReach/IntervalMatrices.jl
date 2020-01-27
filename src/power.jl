import Base: ^

function ^(M::IntervalMatrix{T}, k::Integer) where {T}
    n = checksquare(M)

    # create symbolic matrix
    Ms = Matrix{SymEngine.Basic}(undef, n, n)
    dict = Dict{SymEngine.Basic, Interval{T}}()
    dict2 = Dict{SymEngine.Basic, T}()
    for i in 1:n, j in 1:n
        var = SymEngine.symbols("M_$(i)_$(j)")
        Ms[i, j] = var
        dict[var] = M[i, j]
        dict2[var] = mid(M[i, j])
    end

    # compute symbolic matrix power
    Msᵏ = Ms^k

    # substitute intervals in symbolic result
    Mᵏ = similar(M)
    for i in 1:n, j in 1:n
        expr = Msᵏ[i, j]
#         for (k, v) in dict
#             expr = SymEngine.subs(expr, k, mid(v))
#         end
        expr = SymEngine.subs(expr, dict2...)
        Mᵏ[i, j] = SymEngine.lambdify(expr)
    end
    return Mᵏ
end

using SymEngine, IntervalMatrices
using MacroTools: postwalk

subidx(i) = join(Char.(0x2080 .+ convert.(UInt16, digits(i)[end:-1:1])))

function symbolic_matrix(n; name="M")
    M = Matrix{SymEngine.Basic}(undef, n, n)
    for i in 1:n, j in 1:n
        M[i, j] = name * subidx(i) * subidx(j)
    end
    return M
end

function subs(ex, imat; name="M")
    n = size(imat, 1) # imat assumed square
    for i in 1:n, j in 1:n
        name_ij = Symbol(name * subidx(i) * subidx(j))
        ex == name_ij && return imat[i, j]
    end
    return ex
end

Msym_eval(Msym, M, k) = [postwalk(ex -> subs(ex, M), convert(Expr, m)) for m in Msym]

function power(M, k)
    @assert size(M, 1) == size(M, 2)
    Msym = symbolic_matrix(size(M, 1))
    X = Msym_eval(Msym^k, M, k)
    return IntervalMatrix(eval.(X))
end

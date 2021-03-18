factor = 0.041666666666666664

function Δ(M, i, j, k, l)
    return (M[i, j] - M[i, k]) - (M[l, j] - M[l, k])
end

s(U, X, i, j, k, l) = s(U, X, (i, j, k, l))
function s(U, X, tetrad)
    return Δ(U, tetrad...) * Δ(X, tetrad...)
end

"""
Computes the summand
"""
scombs(X, U, tetrad) = scombs(X, U, tetrad...)
function scombs(X, U, i, j, k, l)
    return factor * (
        Δ(X, i, j, k, l) * U[i,j] +
        Δ(X, i, k, j, l) * (- U[i,j]) +
        Δ(X, k, j, l, i) * (- U[i,j]) +
        Δ(X, l, k, j, i) * U[i,j] +
        Δ(X, k, l, j, i) * U[i,j] +
        Δ(X, l, j, k, i) * (- U[i,j]) +
        Δ(X, i, j, l, k) * U[i,j] +
        Δ(X, i, l, j, k) * (- U[i,j]) 
    )
end

@everywhere Δfe(M, i, j, k, l) = @inbounds M[i, j] - M[i, k] - M[l, j] + M[l, k]

@inline function computeU(X, U, tetrads, peel)

    Nσ = length(tetrads)
    u = 0.
    for part in partition(tetrads, peel)
        
        u += @sync @distributed (+) for t in collect(part)
            Δfe(X, t...) * Δfe(U, t...)
        end
    end

    return u / Nσ

end

@inline function computeU(X, U, tetrads)

    Nσ = length(tetrads)
    u = 0.
    for t in tetrads
        u += Δfe(X, t...) * Δfe(U, t...)
    end

    return u / Nσ

end
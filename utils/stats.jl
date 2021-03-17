factor = 0.041666666666666664

function Δ(M, i, j, k, l)
    return (M[i, j] - M[i, k]) - (M[l, j] - M[l, k])
end

s(U, X, i, j, k, l) = s(U, X, (i, j, k, l))
function s(U, X, tetrad::Tuple{Int64,4})
    return Δ(U, tetrad...) * Δ(X, tetrad...)
end

"""
Computes the summand
"""
scombs(X, U, tetrad) = comb(X, U, tetrad...)
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

"""
Differencing out FE of tetrad of A, using tetrads t
"""
Δfe(M::Matrix) = t -> Δfe(M, t)
Δfe(M::Matrix, t) = Δfe(M::Matrix, t...)
function Δfe(M::Matrix, i, j, k, l)
    first = M[i, j] - M[i, k]
    second = M[l, j] - M[l, k]

    return first - second
end

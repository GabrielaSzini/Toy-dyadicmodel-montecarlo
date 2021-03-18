factor = 0.041666666666666664

function Δ(M, i, j, k, l)
    @inbounds (M[i, j] - M[i, k]) - (M[l, j] - M[l, k])
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
    summand = @inbounds (
        Δ(X, i, j, k, l) * U[i,j] +
        Δ(X, i, k, j, l) * (- U[i,j]) +
        Δ(X, k, j, l, i) * (- U[i,j]) +
        Δ(X, l, k, j, i) * U[i,j] +
        Δ(X, k, l, j, i) * U[i,j] +
        Δ(X, l, j, k, i) * (- U[i,j]) +
        Δ(X, i, j, l, k) * U[i,j] +
        Δ(X, i, l, j, k) * (- U[i,j]) 
    )

    return summand / 24
end

"""
Computes the summand for the second estimator
"""
scombsinefficient(X, U, tetrad) = scombsinefficient(X, U, tetrad...)
function scombsinefficient(X, U, i, j, k, l)
    summands = 0.
    @inbounds for a in (i,j,k,l), b in (i,j,k,l)
        (a - b) == 0 && continue
        for c in (i,j,k,l), d in (i,j,k,l)
            (c - a) * (c - b) == 0 && continue
            (d - c) * (d - b) * (d - a) == 0 && continue
            summands += s(U, X, a, b, c, d) 
        end
    end

    return summands/24

end

Δfe(M, i, j, k, l) = @inbounds M[i, j] - M[i, k] - M[l, j] + M[l, k]

@inline function computeU(
    X::Matrix, U::Matrix, N::Int64)

    Nσ = factorial(N, N - 4)

    u = 0.

    @inbounds for i in 1:N, j in 1:N
        (i - j) == 0 && continue
        for k in 1:N, l in 1:N
            (k - i) * (k - j) == 0 && continue
            (l - k) * (l - j) * (l - i) == 0 && continue
            u += Δfe(X, i, j, k, l) * Δfe(U, i, j, k, l)
        end
    end

    return u / Nσ

end
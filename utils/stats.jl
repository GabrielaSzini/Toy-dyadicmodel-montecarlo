"""
Function that takes double differencing
"""
function Δ(M, i, j, k, l)
    @inbounds M[i, j] - M[i, k] - M[l, j] + M[l, k]
end

"""
Function that computes the score
"""
function s(U, X, i, j, k, l)
    return Δ(U, i, j, k, l) * Δ(X, i, j, k, l)
end

"""
Computes the U-statistic
"""
@inline function computeU(X::Matrix, U::Matrix, N::Int64)

    Nσ = N * (N - 1) * (N - 2) * (N - 3)

    u = 0.

    @inbounds for i in 1:N, j in 1:N
        (i - j) == 0 && continue
        for k in 1:N, l in 1:N
            (k - i) * (k - j) == 0 && continue
            (l - k) * (l - j) * (l - i) == 0 && continue
            u += s(U, X, i, j, k, l)
        end
    end

    return u / Nσ

end

"""
Computes the summand for the estimator of δ₂
"""
function scombspermutations(X, U, i, j, k, l)
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
Computes the summand for the feasible estimator of δ₂
"""
function scombspermutationsfeasible(X, Ũ, i, j, k, l)
    summand = @inbounds (
        Δ(X, i, j, k, l) * Ũ[i,j,k,l] +
        Δ(X, i, k, j, l) * Ũ[i,k,j,l] +
        Δ(X, k, j, l, i) * Ũ[k,j,l,i] +
        Δ(X, l, k, j, i) * Ũ[l,k,j,i] +
        Δ(X, k, l, j, i) * Ũ[k,l,j,i] +
        Δ(X, l, j, k, i) * Ũ[l,j,k,i] +
        Δ(X, i, j, l, k) * Ũ[i,j,l,k] +
        Δ(X, i, l, j, k) * Ũ[i,l,j,k]
    )

    return summand / 24
end

"""
Computes the summand for the estimator of Δ₂
"""
function scombscombinations(X, U, i, j, k, l)
    summand = @inbounds (
        Δ(X, i, j, k, l) * U[i,j] +
        Δ(X, i, k, j, l) * (- U[i,j]) +
        Δ(X, k, j, l, i) * (- U[i,j]) +
        Δ(X, l, k, j, i) * U[i,j] +
        Δ(X, k, l, j, i) * U[i,j] +
        Δ(X, l, j, k, i) * (- U[i,j]) +
        Δ(X, i, j, l, k) * U[i,j] +
        Δ(X, i, l, j, k) * (- U[i,j]) +
        Δ(X, j, i, k, l) * U[j,i] +
        Δ(X, j, k, i, l) * (- U[j,i]) +
        Δ(X, k, i, l, j) * (- U[j,i]) +
        Δ(X, l, k, i, j) * U[j,i] +
        Δ(X, k, l, i, j) * U[j,i] +
        Δ(X, l, i, k, j) * (- U[j,i]) +
        Δ(X, j, i, l, k) * U[j,i] +
        Δ(X, j, l, i, k) * (- U[j,i]) 
    )

    return summand / 24

end

"""
Computes the summand for the feasible estimator of Δ₂
"""
function scombscombinationsfeasible(X, Ũ, i, j, k, l)
    summand = @inbounds (
        Δ(X, i, j, k, l) * Ũ[i, j, k, l] +
        Δ(X, i, k, j, l) * Ũ[i, k, j, l] +
        Δ(X, k, j, l, i) * Ũ[k, j, l, i] +
        Δ(X, l, k, j, i) * Ũ[l, k, j, i] +
        Δ(X, k, l, j, i) * Ũ[k, l, j, i] +
        Δ(X, l, j, k, i) * Ũ[l, j, k, i] +
        Δ(X, i, j, l, k) * Ũ[i, j, l, k] +
        Δ(X, i, l, j, k) * Ũ[i, l, j, k] +
        Δ(X, j, i, k, l) * Ũ[j, i, k, l] +
        Δ(X, j, k, i, l) * Ũ[j, k, i, l] +
        Δ(X, k, i, l, j) * Ũ[k, i, l, j] +
        Δ(X, l, k, i, j) * Ũ[l, k, i, j] +
        Δ(X, k, l, i, j) * Ũ[k, l, i, j] +
        Δ(X, l, i, k, j) * Ũ[l, i, k, j] +
        Δ(X, j, i, l, k) * Ũ[j, i, l, k] +
        Δ(X, j, l, i, k) * Ũ[j, l, i, k] 
    )

    return summand / 24

end

"""
Computes the estimator of δ₂
"""
@inline function computeδ₂(X, U, N)

    Ndyads = N * (N - 1)
    s̄₁ = 0.

    # Fixing a dyad i,j
    ThreadsX.foreach(product(1:N, 1:N)) do (i, j)
        (i - j) == 0 && return

        sdyad₁ = 0.
        Ncombs = binomial(N - 2, 2)

        # Obtaining the conditional expectations over this dyad
        @inbounds for l in 1:N, k in l:N
            (l - i) * (l - j) == 0 && continue
            (l - k) * (k - i) * (k - j) == 0 && continue

            sdyad₁ += scombspermutations(X, U, i, j, k, l)
        end

        s̄₁ += (sdyad₁ / Ncombs)^2
    
    end

    # Taking the average over dyads for both estimators
    δ₂ = s̄₁ / Ndyads
    #Δ₂₂ = s̄₂ / Nyads 

    return δ₂
end

"""
Computes the feasible estimator of δ₂
"""
@inline function computefeasibleδ₂(X, Ũ, N)

    Ndyads = N * (N - 1)
    s̄₁ = 0.

    # Fixing a dyad i,j
    ThreadsX.foreach(product(1:N, 1:N)) do (i, j)
        (i - j) == 0 && return

        sdyad₁ = 0.
        Ncombs = binomial(N - 2, 2)

        # Obtaining the conditional expectations over this dyad
        @inbounds for l in 1:N, k in l:N
            (l - i) * (l - j) == 0 && continue
            (l - k) * (k - i) * (k - j) == 0 && continue

            sdyad₁ += scombspermutationsfeasible(X, Ũ, i, j, k, l)
        end

        s̄₁ += (sdyad₁ / Ncombs)^2
    
    end

    # Taking the average over dyads for both estimators
    δ₂ = s̄₁ / Ndyads
    #Δ₂₂ = s̄₂ / Nyads 

    return δ₂
end


"""
Computes the estimator of Δ₂
"""
@inline function computeΔ₂(X, U, N)

    Ndyadscomb = N * (N - 1) / 2
    s̄₂ = 0.

    # Fixing a combination i,j
    @inbounds for i in 1:N, j in (i+1):N
        (i - j) == 0 && return

        sdyad₂ = 0.
        Ncombs = binomial(N - 2, 2)

        # Obtaining the conditional expectations over this dyad
        @inbounds for l in 1:N, k in l:N
            (l - i) * (l - j) == 0 && continue
            (l - k) * (k - i) * (k - j) == 0 && continue

            sdyad₂ += scombscombinations(X, U, i, j, k, l) 
        end

        s̄₂ += (sdyad₂ / Ncombs)^2
    
    end

    # Taking the average over dyads for both estimators
    Δ₂ = s̄₂ / Ndyadscomb 
    
    return Δ₂
end


"""
Computes the feasible estimator of Δ₂
"""
@inline function computefeasibleΔ₂(X, Ũ, N)

    Ndyadscomb = N * (N - 1) / 2
    s̄₂ = 0.

    # Fixing a combination i,j
    @inbounds for i in 1:N, j in (i+1):N
        (i - j) == 0 && return

        sdyad₂ = 0.
        Ncombs = binomial(N - 2, 2)

        # Obtaining the conditional expectations over this dyad
        @inbounds for l in 1:N, k in l:N
            (l - i) * (l - j) == 0 && continue
            (l - k) * (k - i) * (k - j) == 0 && continue

            sdyad₂ += scombscombinationsfeasible(X, Ũ, i, j, k, l) 
        end

        s̄₂ += (sdyad₂ / Ncombs)^2
    
    end

    # Taking the average over dyads for both estimators
    Δ₂ = s̄₂ / Ndyadscomb 
    
    return Δ₂
end

""" 
Function that computes X̃ or Ỹ for OLS regression
"""
function variablesreg(X::Matrix, N::Int64)
    X̃ = zeros(N * (N-1) * (N-2) * (N-3))

    counter = 1
    # Getting the permutations for constructing the variable
    @inbounds for i in 1:N, j in 1:N # Did not include the threadsX because the order matters
        (i - j) == 0 && continue

        @inbounds for l in 1:N, k in 1:N
            (l - i) * (l - j) == 0 && continue
            (l - k) * (k - i) * (k - j) == 0 && continue

            X̃[counter] = Δ(X, i, j, k, l)
            counter += 1
        end
    end

    return X̃

end

"""
Function that computes the OLS estimator and its variance from standard OLS
"""
function OLSestimator(Y,X)
    estimate = inv(X'*X)*(X'*Y)
    return estimate
end

"""
Function that computes the residuals of the OLS estimation
"""
function Uresiduals(β̂₁, X, Y, N)
    Ũ = zeros(Float64, N, N, N, N)
    @inbounds for i in 1:N, j in 1:N # Did not include the threadsX because the order matters
        (i - j) == 0 && continue

        @inbounds for l in 1:N, k in 1:N
            (l - i) * (l - j) == 0 && continue
            (l - k) * (k - i) * (k - j) == 0 && continue

            Ũ[i,j,k,l] = Δ(Y, i, j, k, l) - β̂₁*Δ(X, i, j, k, l)
        end
    end

    return Ũ
end


"""
Function that computes the variance obtained with U-statistics
"""

function OLSvariance(X, Δ₂, F)
    Mₓₓ = (1/(N * (N-1) * (N-2) * (N-3))) * X' * X
    variance = (1/(N * (N-1))) * inv(Mₓₓ) * F * Δ₂ * inv(Mₓₓ)

    return variance
end

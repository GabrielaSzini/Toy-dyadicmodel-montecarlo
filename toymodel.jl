using DelimitedFiles
using Base.Threads
using Distributed

using Random, Distributions
using Combinatorics
using LinearAlgebra

include("utils/stats.jl")
include("utils/dgp.jl")

function simulation(N)

    E = Set(1:N)

    X, U = datageneration(N)

    Ustat = computeU(X, U, N)
        
    # Obtaining Δ₂
    Ndyads = N * (N - 1)
    s̄₁ = 0.
    s̄₂ = 0.

    # Fixing a dyad i,j
    @inbounds for i in 1:N, j in 1:N
        (i - j) == 0 && continue

        sdyad₁ = 0.
        sdyad₂ = 0.
        Ncombs = binomial(N - 2, 2)

        # Obtaining the conditional expectations over this dyad
        for l in 1:N
            (l - i) * (l - j) == 0 && continue

            for k in l:N
                (k - i) * (k - j) == 0 && continue

                sdyad₁ += scombs(X, U, i, j, k, l)
                sdyad₂ += scombsinefficient(X, U, i, j, k, l) 
            end
        end

        s̄₁ += (sdyad₁ / Ncombs)^2
        s̄₂ += (sdyad₂ / Ncombs)^2
    end

    # Taking the average over dyads for both estimators
    Δ₂₁ = s̄₁ / Ndyads
    Δ₂₂ = s̄₂ / Ndyads 

    return Δ₂₁, Δ₂₂, Ustat 

end

N = 40
sims = 10

@time result = map(simulation, repeat([N], sims))

writedlm("results/out.csv", result, ',')

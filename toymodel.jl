using DelimitedFiles
using Base.Threads

using Random, Distributions
using Combinatorics
using LinearAlgebra

include("utils/stats.jl")
include("utils/dgp.jl")


function simulation(N)

    E = Set(1:N)

    X, U = datageneration(N)

    Ustat = computeU(X, U, N)
        
    # Variance
    Ndyads = N * (N - 1)
    s̄ = 0.

    @inbounds for i in 1:N, j in 1:N
        (i - j) == 0 && continue

        sdyad = 0.
        Ncombs = binomial(N - 2, 2)

        for l in 1:N
            (l - i) * (l - j) == 0 && continue

            for k in l:N
                (k - i) * (k - j) == 0 && continue

                sdyad += scombs(X, U, i, j, k, l)
            end
        end

        s̄ = (sdyad / Ncombs)^2
    end

    Δ₂ = s̄ / Ndyads

    return Δ₂, Ustat 

end

N = 100
sims = 1_000

results = zeros(sims, 2)

@time for i in 1:sims
    print("$i / $sims\r")
    results[i, :] .= simulation(N)
end

writedlm("results/out.csv", results, ',')

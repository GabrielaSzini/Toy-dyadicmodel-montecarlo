using Distributed
using DelimitedFiles

addprocs(8)

@everywhere begin
    using Pkg; Pkg.activate(".")

    using Random, Distributions
    using Combinatorics
    using LinearAlgebra

    include("utils/stats.jl")
    include("utils/dgp.jl")
end

@everywhere function simulation(N)

    E = Set(1:N)

    X, U = datageneration(N)

    tetrads = permutations(1:N, 4)

    Ustat = computeU(X, U, tetrads)
        
        # Variance
    Ndyads = N * (N - 1)
    s̄ = zeros(Ndyads)

    for (l, dyad) in enumerate(permutations(1:N, 2))

        i, j = dyad

            # Combinations of not i, j dyads
        notij = collect(setdiff(E, Set(dyad)))
        othercomb = combinations(notij, 2)

        s̄[l] = mean(scombs(X, U, i, j, k, l) for (k, l) in othercomb)
            
    end

    Δ₂ = mean(s̄.^2)

    return Δ₂, Ustat 

end


N = 100
sims = 1_000

print("Starting parallel simulation...\n")
@time results = pmap(simulation, repeat([N], sims))

print("Not parallel simulation...\n")
@time results = map(simulation, repeat([N], sims))

writedlm("results/out.csv", results, ',')
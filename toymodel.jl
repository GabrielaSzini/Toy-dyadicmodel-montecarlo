using Distributed
using DelimitedFiles
addprocs(4)

@everywhere begin
    using Pkg; Pkg.activate(".")
    using ProgressMeter
    
    using Random, Distributions
    using Combinatorics
    using LinearAlgebra

    include("utils/stats.jl")
    include("utils/dgp.jl")

end


@everywhere function simulation(N, i)

    E = Set(1:N)

    X, U = datageneration(N)

    tetrads = permutations(1:N, 4)

    ΔX, ΔU = Δfe(X), Δfe(U)
    Ustat = mean(@. ΔX(tetrads) * ΔU(tetrads))
    
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

    print("Done with $i\n")

    return Δ₂, Ustat 

end

N = 10
sims = 10

function parallelrun()
    tomap = [(N, i) for _ in 1:sims]
    return pmap(simulation, tomap)
end

@time result = parallelrun()

writedlm("results/out.csv", result, ',')
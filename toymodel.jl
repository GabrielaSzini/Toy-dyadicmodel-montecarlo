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


@everywhere function simulation(N, i)

    E = Set(1:N)

    X, U = datageneration(N)

    tetrads = permutations(1:N, 4)

    ΔX, ΔU = Δfe(X), Δfe(U)
    Ustat = mean(@. ΔX(tetrads) * ΔU(tetrads))
    
    # Variance
    Ndyads = N * (N - 1)
    s̄₁ = zeros(Ndyads)
    s̄₂ = zeros(Ndyads)

    for (l, dyad) in enumerate(permutations(1:N, 2))

        i, j = dyad

        # Combinations of not i, j dyads
        notij = collect(setdiff(E, Set(dyad)))
        othercomb = combinations(notij, 2)

        """
        Here we can propose two estimators for the expectation of s_ij conditional on a dyad i,j
        """     
        s̄₁[l] = mean(scombs(X, U, i, j, k, l) for (k, l) in othercomb)
        s̄₂[l] = mean(scombsinefficient(X, U, i, j, k, l) for (k, l) in othercomb)
    end

    Δ₂₁ = mean(s̄₁.^2)
    Δ₂₂ = mean(s̄₂.^2)

    print("Done with $i\n")

    return  Δ₂₁, Δ₂₂, Ustat 

end

N = 10
sims = 10

function parallelrun()
    tomap = [(N, i) for _ in 1:sims]
    return pmap(simulation, tomap)
end

@time result = parallelrun()

writedlm("results/out.csv", result, ',')
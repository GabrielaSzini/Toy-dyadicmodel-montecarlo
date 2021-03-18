using Distributed
using DelimitedFiles
using Printf
using Plots
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

    return  Δ₂₁, Δ₂₂, Ustat 

end

N = 40
sims = 10

function parallelrun()
    tomap = repeat([N], sims)
    return pmap(simulation, tomap)
end

@time result = parallelrun()

writedlm("results/out.csv", result, ',')

""" 
Some further results to obtain the factors and plots
"""
Δ₂₁ = [x[1] for x in result]
Δ₂₂ = [x[2] for x in result]
Ustat = [x[3] for x in result]

varUstat = var(Ustat)
F₁ = varUstat * (N * (N-1)) ./ Δ₂₁
meanF₁ = mean(F₁)
F̄₁ = varUstat * (N * (N-1)) / mean(Δ₂₁)
@printf("For N = %.0f, and S = %.0f, meanF₁ = %0.5f", N, sims, float(meanF₁))
@printf("For N = %.0f, and S = %.0f, F̄₁ = %0.5f", N, sims, float(F̄₁))

hist₁ = histogram(F₁, bins=:scott, title= "For N = 40, and S = 10", label="", xlabel="F₁",
ylabel="Frequency")
hist₁
savefig("results/histF1_N40_S10.pdf")

F₂ = varUstat * (N * (N-1)) ./ Δ₂₂
meanF₂ = mean(F₂)
F̄₂ = varUstat * (N * (N-1)) / mean(Δ₂₂)
@printf("For N = %.0f, and S = %.0f, meanF₂ = %0.5f", N, sims, float(meanF₂))
@printf("For N = %.0f, and S = %.0f, F̄₂ = %0.5f", N, sims, float(F̄₂))

hist₂ = histogram(F₂, bins=:scott, title= "For N = 40, and S = 10", label="", xlabel="F₂",
ylabel="Frequency")
hist₂
savefig("results/histF2_N40_S10.pdf")
using DelimitedFiles
using Base.Threads
using Printf
using Plots
using Distributed
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

    Ustat = computeU(X, U, N)
        
    # Variance
    Ndyads = N * (N - 1)
    s̄₁ = 0.
    s̄₂ = 0.

    @inbounds for i in 1:N, j in 1:N
        (i - j) == 0 && continue

        sdyad₁ = 0.
        sdyad₂ = 0.
        Ncombs = binomial(N - 2, 2)

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

    Δ₂₁ = s̄₁ / Ndyads
    Δ₂₂ = s̄₂ / Ndyads 

    return Δ₂₁, Δ₂₂, Ustat 

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
using Printf
using Plots
using CSV
using DataFrames

result = DataFrame(CSV.File("results/outN50sims10000.csv", header=0))
""" 
Some further results to obtain the factors and plots
"""
Δ₂₁ = result[:,1]
Δ₂₂ = result[:,2]
Ustat = result[:,3]

varUstat = var(Ustat)
F₁ = varUstat * (N * (N - 1)) ./ Δ₂₁
meanF₁ = mean(F₁)
F̄₁ = varUstat * (N * (N - 1)) / mean(Δ₂₁)
@printf("For N = %.0f, and S = %.0f, meanF₁ = %0.5f", N, sims, float(meanF₁))
@printf("For N = %.0f, and S = %.0f, F̄₁ = %0.5f", N, sims, float(F̄₁))

hist₁ = histogram(F₁, bins=:scott, title="For N = 40, and S = 10", label="", xlabel="F₁",
ylabel="Frequency")
hist₁
savefig("results/histF1_N40_S10.pdf")

F₂ = varUstat * (N * (N - 1)) ./ Δ₂₂
meanF₂ = mean(F₂)
F̄₂ = varUstat * (N * (N - 1)) / mean(Δ₂₂)
@printf("For N = %.0f, and S = %.0f, meanF₂ = %0.5f", N, sims, float(meanF₂))
@printf("For N = %.0f, and S = %.0f, F̄₂ = %0.5f", N, sims, float(F̄₂))

hist₂ = histogram(F₂, bins=:scott, title="For N = 40, and S = 10", label="", xlabel="F₂",
ylabel="Frequency")
hist₂
savefig("results/histF2_N40_S10.pdf")
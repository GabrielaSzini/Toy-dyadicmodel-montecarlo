using Printf
using Plots
using CSV
using DataFrames
using StatsPlots
using Distributions

result = DataFrame(CSV.File("results/outN50sims10000.csv", header=0))
""" 
Some further results to obtain the factors and plots
"""
N=50
sims=10000
β₁=1

Δ₂₁ = result[:,1]
Δ₂₂ = result[:,2]
Ustat = result[:,3]
β̂₁ = result[:,4]
varβ̂₁eff144 = result[:,5]
varβ̂₁eff72 = result[:,6]
varβ̂₁ineff144 = result[:,7]
varβ̂₁ineff72 = result[:,8]

"""
Results for factor analysis
"""    
varUstat = var(Ustat)
F₁ = varUstat * (N * (N - 1)) ./ Δ₂₁
meanF₁ = mean(F₁)
F̄₁ = varUstat * (N * (N - 1)) / mean(Δ₂₁) # look at this one
@printf("For N = %.0f, and S = %.0f, meanF₁ = %0.5f", N, sims, float(meanF₁))
@printf("For N = %.0f, and S = %.0f, F̄₁ = %0.5f", N, sims, float(F̄₁))

hist₁ = histogram(F₁, bins=:scott, title="For N = 50, and S = 10000", label="", xlabel="F₁",
ylabel="Frequency")
hist₁
savefig("results/histograms-F1/histF1_N50_S10000.pdf")

F₂ = varUstat * (N * (N - 1)) ./ Δ₂₂
meanF₂ = mean(F₂)
F̄₂ = varUstat * (N * (N - 1)) / mean(Δ₂₂)
@printf("For N = %.0f, and S = %.0f, meanF₂ = %0.5f", N, sims, float(meanF₂))
@printf("For N = %.0f, and S = %.0f, F̄₂ = %0.5f", N, sims, float(F̄₂))

hist₂ = histogram(F₂, bins=:scott, title="For N = 50, and S = 10000", label="", xlabel="F₂",
ylabel="Frequency")
hist₂
savefig("results/histograms-F2/histF2_N50_S10000.pdf")

"""
Results for OLS estimator
"""

meanβ̂₁ = mean(β̂₁)
@printf("For N = %.0f, and S = %.0f, meanβ̂₁ = %0.5f", N, sims, float(meanβ̂₁))

hist₃ = histogram(β̂₁, bins=:scott, title="For N = 50, and S = 10000", label="", xlabel="β̂₁",
ylabel="Frequency")
hist₃
savefig("results/histograms-beta/histbeta_N50_S10000.pdf")

function myqqplot(obs,F⁰,title)
    nobs=length(obs)
    sort!(obs)
    quantiles⁰ = [quantile(F⁰,i/nobs) for i in 1:nobs]
    # Note that only n-1 points may be plotted, as quantile(F⁰,1) may be inf
    plot(quantiles⁰, obs, seriestype=:scatter, xlabel="Theoretical Quantiles", ylabel = "Sample Quantiles", title=title, label="" )
    plot!(obs,obs,label="")
end
myqqplot(β̂₁,Normal(mean(β̂₁), std(β̂₁)),"For N = 50, and S = 10000")
savefig("results/qqplot-beta/qqplotbeta_N50_S10000.pdf")

varβ̂₁ = var(β̂₁)
@printf("For N = %.0f, and S = %.0f, varβ̂₁ = %0.5f", N, sims, float(varβ̂₁))

meanvarβ̂₁eff144 = mean(varβ̂₁eff144)
@printf("For N = %.0f, and S = %.0f, meanvarβ̂₁eff144 = %0.5f", N, sims, float(meanvarβ̂₁eff144))
sizeβ̂₁eff144 = sum(abs.((β̂₁ .- β₁)./(sqrt.(varβ̂₁eff144))) .> 1.96)/sims 
@printf("For N = %.0f, and S = %.0f, sizeβ̂₁eff144 = %0.5f", N, sims, float(sizeβ̂₁eff144))

meanvarβ̂₁eff72 = mean(varβ̂₁eff72)
@printf("For N = %.0f, and S = %.0f, meanvarβ̂₁eff72 = %0.5f", N, sims, float(meanvarβ̂₁eff72))
sizeβ̂₁eff72 = sum(abs.((β̂₁ .- β₁)./(sqrt.(varβ̂₁eff72))) .> 1.96)/sims 
@printf("For N = %.0f, and S = %.0f, sizeβ̂₁eff72 = %0.5f", N, sims, float(sizeβ̂₁eff72))

meanvarβ̂₁ineff144 = mean(varβ̂₁ineff144)
@printf("For N = %.0f, and S = %.0f, meanvarβ̂₁ineff144 = %0.5f", N, sims, float(meanvarβ̂₁ineff144))
sizeβ̂₁ineff144 = sum(abs.((β̂₁ .- β₁)./(sqrt.(varβ̂₁ineff144))) .> 1.96)/sims 
@printf("For N = %.0f, and S = %.0f, sizeβ̂₁ineff144 = %0.5f", N, sims, float(sizeβ̂₁ineff144))

meanvarβ̂₁ineff72 = mean(varβ̂₁ineff72)
@printf("For N = %.0f, and S = %.0f, meanvarβ̂₁ineff72 = %0.5f", N, sims, float(meanvarβ̂₁ineff72))
sizeβ̂₁ineff72  = sum(abs.((β̂₁ .- β₁)./(sqrt.(varβ̂₁ineff72))) .> 1.96)/sims 
@printf("For N = %.0f, and S = %.0f, sizeβ̂₁ineff72 = %0.5f", N, sims, float(sizeβ̂₁ineff72))

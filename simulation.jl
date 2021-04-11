using Distributed

addprocs(8; exeflags="--project")

@everywhere begin
    using DelimitedFiles
    using ProgressMeter
    using ThreadsX
    using Base.Iterators

    using Random, Distributions
    using LinearAlgebra


    include("utils/stats.jl")
    include("utils/dgp.jl")

end

@everywhere function simulation(N, β₁, design)

    Y, X, U = dgpolsdesign2(N, β₁)

    # Computing U-statistic, δ₂ and Δ₂
    Ustat = computeU(X, U, N)
    δ₂ = computeδ₂(X, U, N)
    Δ₂ = computeΔ₂(X, U, N)

    # Computing OLS estimations
    X̃ = variablesreg(X, N)
    Ỹ = variablesreg(Y, N)
    β̂₁ = OLSestimator(Ỹ,X̃)
    varβ̂₁144 = OLSvariance(X̃, δ₂, 144)
    varβ̂₁72 = OLSvariance(X̃, Δ₂, 72)

    return δ₂, Δ₂, Ustat, β̂₁, varβ̂₁144, varβ̂₁72

end

β₁ = 1
N = 50
sims = 10000
design = 2

result = @time @showprogress pmap(1:sims) do sim
    simulation(N, β₁, design)
end

writedlm("results/outN$(N)sims$(sims)_design$(design).csv", result, ',')
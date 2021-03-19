using Distributed

addprocs(7; exeflags="--project")

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

@everywhere function simulation(N)

    X, U = datageneration(N)

    Ustat = computeU(X, U, N)
    Δ₂₁, Δ₂₂ = computeΔ₂(X, U, N)

    return Δ₂₁, Δ₂₂, Ustat 

end

N = 100
sims = 1000

result = @time @showprogress pmap(1:sims) do sim
    simulation(N)
end

writedlm("results/out.csv", result, ',')
using ThreadsX
using Base.Iterators

using Random, Distributions
using LinearAlgebra

using DelimitedFiles

include("utils/stats.jl")
include("utils/dgp.jl")

function simulation(N)

    E = Set(1:N)

    X, U = datageneration(N)

    @time Ustat = computeU(X, U, N)
    @time Δ₂₁, Δ₂₂ = computeΔ₂(X, U, N)

    return Δ₂₁, Δ₂₂, Ustat 

end

N = 100
sims = 10

@time result = map(simulation, repeat([N], sims))

writedlm("results/out.csv", result, ',')

using Distributed
addprocs(6; exeflags="--project")

@everywhere begin
    using Distributions
    using ForwardDiff
    using StatsFuns
    #using GLM
    using StatFiles
    using DataFrames
    #using MATLAB
    using DelimitedFiles
    using ProgressMeter
    #using ThreadsX
    using Base.Iterators
    using Random
    using LinearAlgebra
    using Optim, NLSolversBase
    #using Base.Threads
    using FileIO, JLD2

    include("utils/gdp_jochmans.jl")
    include("utils/stats.jl")
    include("utils/stats_gdp_jochmans.jl")

end


@everywhere function simulation(N, θ₀, Cₙ, dimension_reg)

    #DGP of Jochmans
    Y, X = gdp_jochmans(N, θ₀, Cₙ)

    index_i_comb, index_j_comb, ỹ_comb, b_comb, c_comb, nondiag_comb, X̃_comb =  dataquantilecomb_gdp_jochmans(Y, X)

    X̃_comb_cond = @view X̃_comb[c_comb.==1]
    ỹ_comb_cond = @view ỹ_comb[c_comb.==1]
    nondiag_comb_cond = @view nondiag_comb[c_comb.==1]

    β₀_comb = zeros(1,1)
    β₀_comb = initialize!(β₀_comb, X̃_comb_cond, ỹ_comb_cond)
    nll_comb = make_closures(X̃_comb_cond, ỹ_comb_cond)
    result_comb = optimize(nll_comb, β₀_comb, LBFGS(), autodiff=:forward)
    β̂_comb = Optim.minimizer(result_comb)

    # point estimates with permutations
    index_i_perm, index_j_perm, ỹ_perm, b_perm, c_perm, nondiag_perm, X̃_perm = dataquantileperm_gdp_jochmans(Y, X)

    X̃_perm_cond = @view X̃_perm[c_perm.==1]
    ỹ_perm_cond = @view ỹ_perm[c_perm.==1]

    β₀_perm = zeros(1,1)
    β₀_perm = initialize!(β₀_perm, X̃_perm_cond, ỹ_perm_cond)
    nll_perm = make_closures(X̃_perm_cond, ỹ_perm_cond)
    result_perm = optimize(nll_perm, β₀_perm, LBFGS(), autodiff=:forward)
    β̂_perm = Optim.minimizer(result_perm)

    #standard errors
    se₁, se₂, se₃, se₄ = standarderrors_gdp_jochmans_new(Y, X, β̂_comb, X̃_comb_cond, nondiag_comb_cond)

    Y[diagind(Y)] .= 0.0
    perc_links = sum(Y)/(N*(N-1))
    perc_quadruples = sum(c_comb)/length(c_comb)

    return β̂_comb, β̂_perm, se₁, se₂, se₃, se₄, perc_quadruples, perc_links

end

########### N = 30, DGP 1 ##########
#Defining number of nodes
N=30
#Defining the parameter
θ₀ = 1
#Defining sparsity
Cₙ = 0
#Cₙ = log(log(N))
#Cₙ = (log(N))^(1/2)
#Cₙ = log(N)
#Defining dimension_reg
dimension_reg = 1
#Defining number of simulations
sims = 1000

result1_30 = @time @showprogress pmap(1:sims) do sim
    simulation(N, θ₀, Cₙ, dimension_reg)
end

println("Finalized results DGP 1, N=30 - Jochmans")

save("results_simulations/results_jochmans_DGP1_N30.jld2","result1_30",result1_30)

########### N = 30, DGP 2 ##########
#Defining number of nodes
N=30
#Defining the parameter
θ₀ = 1
#Defining sparsity
#Cₙ = 0
Cₙ = log(log(N))
#Cₙ = (log(N))^(1/2)
#Cₙ = log(N)
#Defining dimension_reg
dimension_reg = 1
#Defining number of simulations
sims = 1000

result2_30 = @time @showprogress pmap(1:sims) do sim
    simulation(N, θ₀, Cₙ, dimension_reg)
end

println("Finalized results DGP 2, N=30 - Jochmans")

save("results_simulations/results_jochmans_DGP2_N30.jld2","result2_30",result2_30)

########### N = 30, DGP 3 ##########
#Defining number of nodes
N=30
#Defining the parameter
θ₀ = 1
#Defining sparsity
#Cₙ = 0
#Cₙ = log(log(N))
Cₙ = (log(N))^(1/2)
#Cₙ = log(N)
#Defining dimension_reg
dimension_reg = 1
#Defining number of simulations
sims = 1000

result3_30 = @time @showprogress pmap(1:sims) do sim
    simulation(N, θ₀, Cₙ, dimension_reg)
end

println("Finalized results DGP 3, N=30 - Jochmans")

save("results_simulations/results_jochmans_DGP3_N30.jld2","result3_30",result3_30)

########### N = 30, DGP 4 ##########
#Defining number of nodes
N=30
#Defining the parameter
θ₀ = 1
#Defining sparsity
#Cₙ = 0
#Cₙ = log(log(N))
#Cₙ = (log(N))^(1/2)
Cₙ = log(N)
#Defining dimension_reg
dimension_reg = 1
#Defining number of simulations
sims = 1000

result4_30 = @time @showprogress pmap(1:sims) do sim
    simulation(N, θ₀, Cₙ, dimension_reg)
end

println("Finalized results DGP 4, N=30 - Jochmans")

save("results_simulations/results_jochmans_DGP4_N30.jld2","result4_30",result4_30)

########### N = 45, DGP 1 ##########
#Defining number of nodes
N=45
#Defining the parameter
θ₀ = 1
#Defining sparsity
Cₙ = 0
#Cₙ = log(log(N))
#Cₙ = (log(N))^(1/2)
#Cₙ = log(N)
#Defining dimension_reg
dimension_reg = 1
#Defining number of simulations
sims = 1000

result1_45 = @time @showprogress pmap(1:sims) do sim
    simulation(N, θ₀, Cₙ, dimension_reg)
end

println("Finalized results DGP 1, N=45 - Jochmans")

save("results_simulations/results_jochmans_DGP1_N45.jld2","result1_45",result1_45)

########### N = 45, DGP 2 ##########
#Defining number of nodes
N=45
#Defining the parameter
θ₀ = 1
#Defining sparsity
#Cₙ = 0
Cₙ = log(log(N))
#Cₙ = (log(N))^(1/2)
#Cₙ = log(N)
#Defining dimension_reg
dimension_reg = 1
#Defining number of simulations
sims = 1000

result2_45 = @time @showprogress pmap(1:sims) do sim
    simulation(N, θ₀, Cₙ, dimension_reg)
end

println("Finalized results DGP 2, N=45 - Jochmans")

save("results_simulations/results_jochmans_DGP2_N45.jld2","result2_45",result2_45)

########### N = 45, DGP 3 ##########
#Defining number of nodes
N=45
#Defining the parameter
θ₀ = 1
#Defining sparsity
#Cₙ = 0
#Cₙ = log(log(N))
Cₙ = (log(N))^(1/2)
#Cₙ = log(N)
#Defining dimension_reg
dimension_reg = 1
#Defining number of simulations
sims = 1000

result3_45 = @time @showprogress pmap(1:sims) do sim
    simulation(N, θ₀, Cₙ, dimension_reg)
end

println("Finalized results DGP 3, N=45 - Jochmans")

save("results_simulations/results_jochmans_DGP3_N45.jld2","result3_45",result3_45)

########### N = 30, DGP 4 ##########
#Defining number of nodes
N=45
#Defining the parameter
θ₀ = 1
#Defining sparsity
#Cₙ = 0
#Cₙ = log(log(N))
#Cₙ = (log(N))^(1/2)
Cₙ = log(N)
#Defining dimension_reg
dimension_reg = 1
#Defining number of simulations
sims = 1000

result4_45 = @time @showprogress pmap(1:sims) do sim
    simulation(N, θ₀, Cₙ, dimension_reg)
end

println("Finalized results DGP 4, N=45 - Jochmans")

save("results_simulations/results_jochmans_DGP4_N45.jld2","result4_45",result4_45)
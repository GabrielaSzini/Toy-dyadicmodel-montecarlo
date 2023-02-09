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

    include("utils/gdp_dzemski.jl")
    include("utils/stats.jl")
    #note that we can use the same codes as in Jochmans!
    include("utils/stats_gdp_jochmans.jl")

    function gdp_jochmans_fe(N, θ₀, Cₙ)

        #Generating regressor xᵢⱼ
        uᵢ = rand(Beta(2,2),N) .- 1/2
        uᵢ = reshape(uᵢ,N,1)
        xᵢ = repeat(uᵢ, 1, N)
        uⱼ = uᵢ'
        xⱼ = repeat(uⱼ,N,1)
        xᵢⱼ = - abs.(xᵢ .- xⱼ)
    
        #Generate fixed effects
        αᵢ = - (N .- reshape([1:1:N;],N,1))./(N-1).*Cₙ
        αᵢ = repeat(αᵢ, 1,N)
        γⱼ = αᵢ'
    
        #Generate error term
        ϵᵢⱼ = rand(Logistic(0,1),N,N)
    
        #Generate dependent variable
        yᵢⱼ = Int.(xᵢⱼ.*θ₀ .+ αᵢ .+ γⱼ .- ϵᵢⱼ .≥ 0)
    
        return αᵢ
    
    end

    function gdp_dzemski_fe(N, θ₀, Cₙ)

        #Generating regressor xᵢⱼ
        xᵢ = ones(N) - 2 .* isodd.(reshape([1:1:N;],N,1))
        xᵢ = repeat(xᵢ, 1, N)
        xⱼ = xᵢ'
        xᵢⱼ = xᵢ .* xⱼ
    
        #Generate fixed effects
        αᵢ = - (N .- reshape([1:1:N;],N,1))./(N-1).*Cₙ
        αᵢ = repeat(αᵢ, 1,N)
        γⱼ = αᵢ'
    
        #Generate error term
        ϵᵢⱼ = rand(Logistic(0,1),N,N)
    
        #Generate dependent variable
        yᵢⱼ = Int.(xᵢⱼ.*θ₀ .+ αᵢ .+ γⱼ .- ϵᵢⱼ .≥ 0)
    
        return αᵢ
    
    end

    function gdp_yan(N, θ₀₁, θ₀₂, Cₙ)

        #Generating regressor xᵢⱼ
        uᵢ₁ = rand(Beta(2,2),N) 
        uᵢ₁ = reshape(uᵢ₁,N,1)
        xᵢ₁ = repeat(uᵢ₁, 1, N)
        uⱼ₁ = uᵢ₁'
        xⱼ₁ = repeat(uⱼ₁,N,1)
        xᵢⱼ₁ = - abs.(xᵢ₁ .- xⱼ₁)
    
        uᵢ₂ = rand(Beta(2,2),N) 
        uᵢ₂ = reshape(uᵢ₂,N,1)
        xᵢ₂ = repeat(uᵢ₂, 1, N)
        uⱼ₂ = uᵢ₂'
        xⱼ₂ = repeat(uⱼ₂,N,1)
        xᵢⱼ₂ = - abs.(xᵢ₂ .- xⱼ₂)
    
    
        #Generate fixed effects
        αᵢ = - (N .- reshape([1:1:N;],N,1))./(N-1).*Cₙ
        αᵢ = repeat(αᵢ, 1,N)
        γⱼ = αᵢ'
    
        #Generate error term
        ϵᵢⱼ = rand(Logistic(0,1),N,N)
    
        #Generate dependent variable
        yᵢⱼ = Int.(xᵢⱼ₁.*θ₀₁ .+ xᵢⱼ₂.*θ₀₂ .+ αᵢ .+ γⱼ .- ϵᵢⱼ .≥ 0)
    
        return αᵢ
    
    end

    N = 10000
    θ₀ = 1
    θ₀₁ = 1
    θ₀₂ = 1.5

    Cₙ = 0
    ß = gdp_jochmans_fe(N, θ₀, Cₙ)
    fe_dzemski_cn1 = gdp_dzemski_fe(N, θ₀, Cₙ)
    fe_yan_cn1 = gdp_yan(N, θ₀₁, θ₀₂, Cₙ)

    Cₙ = log(log(N))
    fe_jochmans_cn2 = gdp_jochmans_fe(N, θ₀, Cₙ)
    fe_dzemski_cn2 = gdp_dzemski_fe(N, θ₀, Cₙ)
    fe_yan_cn2 = gdp_yan(N, θ₀₁, θ₀₂, Cₙ)

    Cₙ = (log(N))^(1/2)
    fe_jochmans_cn3 = gdp_jochmans_fe(N, θ₀, Cₙ)
    fe_dzemski_cn3 = gdp_dzemski_fe(N, θ₀, Cₙ)
    fe_yan_cn3 = gdp_yan(N, θ₀₁, θ₀₂, Cₙ)

    Cₙ = log(N)
    fe_jochmans_cn4 = gdp_jochmans_fe(N, θ₀, Cₙ)
    fe_dzemski_cn4 = gdp_dzemski_fe(N, θ₀, Cₙ)
    fe_yan_cn4 = gdp_yan(N, θ₀₁, θ₀₂, Cₙ)
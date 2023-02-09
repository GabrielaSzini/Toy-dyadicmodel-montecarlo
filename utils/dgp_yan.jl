function dgp_yan(N, θ₀₁, θ₀₂, Cₙ)

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

    return yᵢⱼ, xᵢⱼ₁, xᵢⱼ₂

end
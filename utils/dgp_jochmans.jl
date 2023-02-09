function dgp_jochmans(N, θ₀, Cₙ)

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

    return yᵢⱼ, xᵢⱼ

end
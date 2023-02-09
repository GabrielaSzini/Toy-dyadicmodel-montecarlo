function dgp_dzemski(N, θ₀, Cₙ)

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

    return yᵢⱼ, xᵢⱼ

end
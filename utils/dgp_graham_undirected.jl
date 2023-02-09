function dgp_graham_undirected(N, θ₀, αₗ, αₕ, λ₀, λ₁)

    #Generating regressor xᵢⱼ
    xᵢ = rand(Normal(0,1),N)
    xᵢ[xᵢ.<=0] .= -1
    xᵢ[xᵢ.> 0] .= 1
    xᵢ = reshape(xᵢ,N,1)
    xᵢ = repeat(xᵢ, 1, N)
    xⱼ = xᵢ'
    xᵢⱼ = xᵢ .* xⱼ

    #Generate fixed effects
    condₗ = xᵢ .== -1
    condₕ = xᵢ .== 1
    vᵢ = rand(Beta(λ₀,λ₁),N) .- λ₀/(λ₀ + λ₁)
    vᵢ = reshape(vᵢ,N,1)
    vᵢ = repeat(vᵢ, 1, N)
    aᵢ = αₗ .* condₗ .+ αₕ .* condₕ .+ vᵢ
    aⱼ = aᵢ'

    #Generate error term
    ϵᵢⱼ_directed = rand(Logistic(0,1),N,N)
    ϵᵢⱼ = Symmetric(ϵᵢⱼ_directed, :L)


    #Generate dependent variable
    yᵢⱼ = Int.(xᵢⱼ.*θ₀ .+ aᵢ .+ aⱼ .- ϵᵢⱼ .≥ 0)

    return yᵢⱼ, xᵢⱼ, aᵢ, aⱼ

end

# if we want the network to be undirected, then the error term should be symmetric!
"""
Generates the X and U
"""
function datageneration(N)

    Uij = rand(Normal(), (N, N))
    Uij[diagind(Uij)] .= 0.

    # drawing the FE likewise Jochmans (2016)
    Ai = rand(Beta(2, 2), (N, 1)) .- 0.5
    Xij = @. - abs(Ai - Ai')

    return Xij, Uij

end 

"""
DGP for OLS estimator: first design (fixed effects uncorrelated with explanatory variables and continuous X)
"""

function dgpolsfirst(N, β₁)

    Uij = rand(Normal(), (N, N))
    Uij[diagind(Uij)] .= 0.

    # drawing the FE for the explanatory variables likewise Jochmans (2016)
    Ai = rand(Beta(2, 2), (N, 1)) .- 0.5
    Xij = @. - abs(Ai - Ai')
    Xij[diagind(Xij)] .= 0.

    # drawing the uncorrelated FE for the regression
    θᵢ = rand(Normal(), N)
    ξⱼ = rand(Normal(), N)

    θ = repeat(θᵢ, 1, N) # in each row we have the repeated FE for individual i
    ξ = repeat(ξⱼ, 1, N)' # in each column we have the repeated FE for individual j

    Yij = β₁ * Xij + θ + ξ + Uij
    Yij[diagind(Yij)] .= 0.

    return Yij, Xij, Uij
end

"""
DGP for OLS estimator: second design (fixed effects correlated with explanatory variables and continuous X)
""" 

function dgpolssecond(N, β₁)

    Uij = rand(Normal(), (N, N))
    Uij[diagind(Uij)] .= 0.

    # drawing the FE for the regression that will be correlated to the explanatory variables
    θᵢ = rand(Normal(), N)
    ξⱼ = rand(Normal(), N)

    θ = repeat(θᵢ, 1, N) # in each row we have the repeated FE for individual i
    ξ = repeat(ξⱼ, 1, N)' # in each column we have the repeated FE for individual j

    # drawing the FE for the explanatory variables likewise Jochmans (2016), but including Aj and correlation with θ and ξ
    Ai = rand(Beta(2, 2), (N, 1)) .- 0.5
    Aj = rand(Beta(2, 2), (N, 1)) .- 0.5
    Xij = @. - abs(Ai - Aj')
    Xij = Xij + θ + ξ
    Xij[diagind(Xij)] .= 0.   

    Yij = β₁ * Xij + θ + ξ + Uij
    Yij[diagind(Yij)] .= 0.

    return Yij, Xij, Uij
end

"""
DGP for OLS estimator: third design (fixed effects uncorrelated with explanatory variables and discrete X)
"""
function dgpolsthird(N, β₁)

    Uij = rand(Normal(), (N, N))
    Uij[diagind(Uij)] .= 0.

    # drawing the FE for the regression that will be correlated to the explanatory variables
    θᵢ = rand(Normal(), N)
    ξⱼ = rand(Normal(), N)

    θ = repeat(θᵢ, 1, N) # in each row we have the repeated FE for individual i
    ξ = repeat(ξⱼ, 1, N)' # in each column we have the repeated FE for individual j

    # drawing the FE for the explanatory variables likewise Jochmans (2016), but including Aj and correlation with θ and ξ
    Ai = rand(Beta(2, 2), (N, 1)) .- 0.5
    Aj = rand(Beta(2, 2), (N, 1)) .- 0.5
    Xij = @. (Ai - Aj')
    Xij = (Xij) .> 0
    Xij = convert(Array{Float64}, Xij)
    Xij[diagind(Xij)] .= 0.   

    Yij = β₁ * Xij + θ + ξ + Uij
    Yij[diagind(Yij)] .= 0.

    return Yij, Xij, Uij
end

"""
DGP for OLS estimator: fourth design (fixed effects correlated with explanatory variables and discrete X)
"""
function dgpolsfourth(N, β₁)

    Uij = rand(Normal(), (N, N))
    Uij[diagind(Uij)] .= 0.

    # drawing the FE for the regression that will be correlated to the explanatory variables
    θᵢ = rand(Normal(), N)
    ξⱼ = rand(Normal(), N)

    θ = repeat(θᵢ, 1, N) # in each row we have the repeated FE for individual i
    ξ = repeat(ξⱼ, 1, N)' # in each column we have the repeated FE for individual j

    # drawing the FE for the explanatory variables likewise Jochmans (2016), but including Aj and correlation with θ and ξ
    Ai = rand(Beta(2, 2), (N, 1)) .- 0.5
    Aj = rand(Beta(2, 2), (N, 1)) .- 0.5
    Xij = @. (Ai - Aj')
    Xij = (Xij + θ + ξ) .> 0
    Xij = convert(Array{Float64}, Xij)
    Xij[diagind(Xij)] .= 0.   

    Yij = β₁ * Xij + θ + ξ + Uij
    Yij[diagind(Yij)] .= 0.

    return Yij, Xij, Uij
end

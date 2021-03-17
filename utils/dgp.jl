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
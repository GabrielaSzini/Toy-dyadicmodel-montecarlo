include("utils/dgp_graham.jl")
include("utils/dgp_graham_undirected.jl")
using Random
using Combinatorics
using Distributions
using Base.Iterators
using JuMP
using Clp
using HiGHS
using LinearAlgebra

grid_points = 5
N = 100
αₗ = -0.5
αₕ = -0.5
λ₀ = 1
λ₁ = 1
θ₀ = 1
subnet_size = 4

# Take simulation
Y, X, Aᵢ, Aⱼ = dgp_graham_undirected(N, θ₀, αₗ, αₕ, λ₀, λ₁)

# Create grids for aᵢ
minimum_a = αₗ - (λ₀/(λ₁ + λ₀))
maximum_a = αₕ + (λ₀/(λ₁ + λ₀))
a_grid = LinRange(minimum_a,maximum_a,grid_points)
total_a_grid = collect(Iterators.product(a_grid,a_grid,a_grid,a_grid))
total_a_grid = reshape(total_a_grid,grid_points^subnet_size,1)
a_1 = [ x[1] for x in total_a_grid ]
a_2 = [ x[2] for x in total_a_grid ]
a_3 = [ x[3] for x in total_a_grid ]
a_4 = [ x[4] for x in total_a_grid ]
# get combination of fixed effects for each combination of total_grid_a => rows are combinations, 625 of those => matrix 625x6
a_i_matrix= hcat(a_1, a_1, a_1, a_2, a_2, a_3)
a_j_matrix= hcat(a_2, a_3, a_4, a_3, a_4, a_4)

# as we are now working with undirected networks, take only the upper triangular matrix of the dgp
# Create combinations of y's : y_12, y_13, y_14, y_23, y_24, y_34, that is, subnetsize*(subnetsize-1)/2
y_possible = (0,1)
total_y_grid = collect(Iterators.product(y_possible,y_possible,y_possible,y_possible,y_possible,y_possible))
total_y_grid = reshape(total_y_grid,Int(2^(subnet_size*(subnet_size-1)/2)),1)
total_y_grid_matrix = reshape(vec(reduce(hcat, getindex.(total_y_grid,i) for i in eachindex(total_y_grid[1]))),64,6)

# Change with respect to main file is that now we do not take random subnetworks, but according to tetrad speficication of Graham
# Subset subnetworks

#nodes = shuffle(1:N)
#s = Int(N/subnet_size)
#subnetworks = collect(partition(nodes,subnet_size))

nodes = collect(1:N)
#using the function combinations tend not to work for high values of N
#combination_nodes =  collect(combinations(nodes,4))

#obtaining quadruples

combination_1 = []
combination_2 = []
combination_3 = []
combination_4 = []

@inbounds for i in 1:N, j in (i+1):N
    (i - j) == 0 && continue
    for k in (j+1):N, l in (k+1):N
        (k - i) * (k - j) == 0 && continue
        (l - k) * (l - j) * (l - i) == 0 && continue

        # got a quadruple, now need to look at the permutations and the Sijk evalueated at those
        combination_nodes = [i,j,k,l]
        permutation_nodes = collect(permutations(combination_nodes))
        include_comb = zeros(length(permutation_nodes))
        for perm in 1:length(permutation_nodes)
            si = permutation_nodes[perm][1]
            sj = permutation_nodes[perm][2]
            sk = permutation_nodes[perm][3]
            sl = permutation_nodes[perm][4]

            Sijkl = Y[si,sj]*Y[sk,sl]*(1-Y[si,sk])*(1-Y[sj,sl]) - (1- Y[si,sj])*(1-Y[sk,sl])*Y[si,sk]*Y[sj,sl]
            if Sijkl == 1 ||  Sijkl == -1
                include_comb[perm] = 1
            else 
                include_comb[perm] = 0
            end
        end

        if sum(include_comb) != 0 
            push!(combination_1, i)
            push!(combination_2, j)
            push!(combination_3, k)
            push!(combination_4, l)
        end
    end
end


# Taking x vector as fixed : x_12, x_13, x_14, x_23, x_24, x_34
fix_x_matrix = Array{NTuple{6,Float64}}(undef, length(combination_1))
for i in 1:length(combination_1)
    ind1 = combination_1[i]
    ind2 = combination_2[i]
    ind3 = combination_3[i]
    ind4 = combination_4[i]
    fix_x_matrix[i] = (X[ind1, ind2], X[ind1, ind3], X[ind1, ind4], X[ind2, ind3], X[ind2, ind4], X[ind3, ind4])
end

fix_y_matrix = Array{NTuple{6,Float64}}(undef, length(combination_1))
for i in 1:length(combination_1)
    ind1 = combination_1[i]
    ind2 = combination_2[i]
    ind3 = combination_3[i]
    ind4 = combination_4[i]
    fix_y_matrix[i] = (Y[ind1, ind2], Y[ind1, ind3],Y[ind1, ind4], Y[ind2, ind3], Y[ind2, ind4], Y[ind3, ind4])
end

###### First set of conditions

#f(y, aᵢ, aⱼ, x, θ) = 1 / (1 + exp(-(x * θ₀ + a_i + a_j))) 

## Marginal effect
# get changes in covariates for all combinations of grid grid_points
x₁ = ones(length(total_a_grid), 6)
x₂ = -1 .*ones(length(total_a_grid), 6)

m = exp.(x₁.*θ₀ .+ a_i_matrix .+ a_j_matrix)./( 1 .+ exp.(x₁.*θ₀ .+ a_i_matrix .+ a_j_matrix)) .- exp.(x₂.*θ₀ .+ a_i_matrix .+ a_j_matrix)./( 1 .+ exp.(x₂.*θ₀ .+ a_i_matrix .+ a_j_matrix))
m = mean(m, dims=2)

result_l_y = zeros(length(total_y_grid), length(combination_1))
result_u_y = zeros(length(total_y_grid), length(combination_1))

# Fixing one subnetwork: one optimization for each subnetwork
for sub in 1:length(combination_1)
    
    fix_x_vector = fix_x_matrix[sub]

    ## F matrix
    f_matrix = zeros(length(total_y_grid), length(total_a_grid))
    for i in 1:length(total_y_grid)
        for j in 1:length(total_a_grid)
            f_matrix[i,j] =  prod(((exp.(fix_x_vector.*θ₀ .+ a_i_matrix[j,:] .+ a_j_matrix[j,:])).^total_y_grid_matrix[i,:])./( 1 .+ exp.(fix_x_vector.*θ₀ .+ a_i_matrix[j,:] .+ a_j_matrix[j,:]))) 
        end
    end

    ##Checking if it makes sense - this is actually set up inside of JuMP
    #u_y = fill(1 , length(total_y_grid))
    #l_y = fill(-1, length(total_y_grid))
    #checking = (u_y' * f_matrix)'

    ###### Second set of conditions
    bmax = 1
    bmin = -1

    ###### Checking if criterion function makes sense
    #checking = sum((u_y - l_y)'*f_matrix)

    ###### Optimization with JuMP

    n = length(total_y_grid)
    bounds = Model(HiGHS.Optimizer)
    set_silent(bounds)

    # create variables to be optimized over
    @variable(bounds, l_y[1:n] >= bmin)
    @variable(bounds, u_y[1:n] <= bmax)
    # create constraints
    @constraint(bounds, f_matrix'l_y .<= m)
    @constraint(bounds, f_matrix'u_y .>= m)
    @constraint(bounds, [i = 1:n], l_y[i] <= u_y[i])
    # create Objective
    @objective(bounds, Min, sum((u_y - l_y)'f_matrix))
    optimize!(bounds)
    objective_value(bounds)
    
    result_l_y[:,sub] = value.(l_y)
    result_u_y[:,sub] = value.(u_y)

end

# find indices of observed y's in total_y_grid
indices = Array{Array{Int64,1},1}(undef, length(fix_y_matrix))

for col in 1:length(fix_y_matrix)
    row = Tuple.(findall(x->x==fix_y_matrix[col], total_y_grid))[1][1]
    indices[col] = [floor(Int64,row),floor(Int64,col)]
end

result_l_y_real = [result_l_y[x,y] for (x,y) in indices]
result_u_y_real = [result_u_y[x,y] for (x,y) in indices]

result_l_y_mean = mean(result_l_y_real)
result_u_y_mean = mean(result_u_y_real)

result_u_y_mean .>= result_l_y_mean

#### Compute true marginal effects
# we take only the lower triangular elements as we have a directed network

Y_lower = tril(Y,-1)
Aᵢ_lower = tril(Aᵢ,-1)
Aⱼ_lower = tril(Aⱼ, -1)

x₁_lower = ones(size(Y)[1], size(Y)[2])
x₁_lower = tril(x₁_lower,-1)
x₂_lower = -1 .*ones(size(Y)[1], size(Y)[2])
x₂_lower = tril(x₂_lower,-1)

individual_marginal = exp.(x₁_lower.*θ₀ .+ Aᵢ_lower .+ Aⱼ_lower)./( 1 .+ exp.(x₁_lower.*θ₀ .+ Aᵢ_lower .+ Aⱼ_lower)) .- exp.(x₂_lower.*θ₀ .+ Aᵢ_lower .+ Aⱼ_lower)./( 1 .+ exp.(x₂_lower.*θ₀ .+ Aᵢ_lower .+ Aⱼ_lower))
true_marginal = sum(individual_marginal)/((N*(N-1))/2)
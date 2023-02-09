using JuMP

"""
Thinking with 5 nodes all possible non directed connections
12 13 14 15 23 24 25 34 35 45
1  1  1  1  2  2  2  3  3  4
2  3  4  5  3  4  5  4  5  5
"""

"""
Function that creates a grids
"""
function agrids(αₗ, αₕ, λ₀, λ₁, grid_points, subnet_size)
    minimum_a = αₗ - (λ₀/(λ₁ + λ₀))
    maximum_a = αₕ + (λ₀/(λ₁ + λ₀))
    a_grid = LinRange(minimum_a,maximum_a,grid_points)
    total_a_grid = collect(Iterators.product(a_grid,a_grid,a_grid,a_grid,a_grid))
    total_a_grid = reshape(total_a_grid,grid_points^subnet_size,1)
    a_1 = [ x[1] for x in total_a_grid ]
    a_2 = [ x[2] for x in total_a_grid ]
    a_3 = [ x[3] for x in total_a_grid ]
    a_4 = [ x[4] for x in total_a_grid ]
    a_5 = [ x[5] for x in total_a_grid ]
    # get combination of fixed effects for each combination of total_grid_a => rows are combinations, 625 of those => matrix 625x6
    a_i_matrix= hcat(a_1, a_1, a_1, a_1, a_2, a_2, a_2, a_3, a_3, a_4)
    a_j_matrix= hcat(a_2, a_3, a_4, a_5, a_3, a_4, a_5, a_4, a_5, a_5)

    return a_i_matrix, a_j_matrix, total_a_grid

end

"""
Function that creates y grids
"""
function ygrids(subnet_size)
    y_possible = (0,1)
    total_y_grid = collect(Iterators.product(y_possible,y_possible,y_possible,y_possible,y_possible,y_possible,y_possible,y_possible,y_possible,y_possible))
    total_y_grid = reshape(total_y_grid,Int(2^(subnet_size*(subnet_size-1)/2)),1)
    total_y_grid_matrix = reshape(vec(reduce(hcat, getindex.(total_y_grid,i) for i in eachindex(total_y_grid[1]))),1024,10)

    return total_y_grid, total_y_grid_matrix

end

"""
Function that creates matrices of Y and X for random subnetworks 
"""
function fix_yx(subnetworks, Y, X, s)
    # Taking x vector as fixed : x_12, x_13, x_14, x_23, x_24, x_34
    fix_x_matrix = Array{NTuple{10,Float64}}(undef, s)
    fix_y_matrix = Array{NTuple{10,Float64}}(undef, s)
    for i in 1:Int(s)
        ind1 = subnetworks[i][1]
        ind2 = subnetworks[i][2]
        ind3 = subnetworks[i][3]
        ind4 = subnetworks[i][4]
        ind5 = subnetworks[i][5]
        fix_x_matrix[i] = (X[ind1, ind2], X[ind1, ind3], X[ind1, ind4], X[ind1, ind5], X[ind2, ind3], X[ind2, ind4], X[ind2, ind5], X[ind3, ind4], X[ind3, ind5], X[ind4, ind5])
        fix_y_matrix[i] = (Y[ind1, ind2], Y[ind1, ind3], Y[ind1, ind4], Y[ind1, ind5], Y[ind2, ind3], Y[ind2, ind4], Y[ind2, ind5], Y[ind3, ind4], Y[ind3, ind5], Y[ind4, ind5])
    end

    return fix_x_matrix, fix_y_matrix

end

"""
Function that computes the value of m for the first condition (given an x₁ and x₂)
"""
function mvalue(x₁_value, x₂_value, total_a_grid, θ₀, a_i_matrix, a_j_matrix)
    x₁ = x₁_value.*ones(length(total_a_grid), 10)
    x₂ = x₂_value.*ones(length(total_a_grid), 10)

    m = exp.(x₁.*θ₀ .+ a_i_matrix .+ a_j_matrix)./( 1 .+ exp.(x₁.*θ₀ .+ a_i_matrix .+ a_j_matrix)) .- exp.(x₂.*θ₀ .+ a_i_matrix .+ a_j_matrix)./( 1 .+ exp.(x₂.*θ₀ .+ a_i_matrix .+ a_j_matrix))
    m = mean(m, dims=2)

    return m

end

"""
Function that makes optimization for random subnetworks
"""
function mainopt(total_y_grid, subnetworks, s, fix_x_matrix, total_a_grid, θ₀, a_i_matrix, a_j_matrix, total_y_grid_matrix, bmax, bmin, m)

    #result_l_y = zeros(length(total_y_grid), length(combinations_final))
    result_l_y = Array{Union{Float64,Missing}}(missing,length(total_y_grid),s) 
    #result_u_y = zeros(length(total_y_grid), length(combinations_final))
    result_u_y = Array{Union{Float64,Missing}}(missing,length(total_y_grid),s) 

    # Fixing one subnetwork: one optimization for each subnetwork
    for sub in 1:s
    
        fix_x_vector = fix_x_matrix[sub]

        ## F matrix
        #f_matrix = zeros(length(total_y_grid), length(total_a_grid))
        f_matrix = Array{Union{Float64,Missing}}(missing,length(total_y_grid),length(total_a_grid)) 

        for i in 1:length(total_y_grid)
            for j in 1:length(total_a_grid)
                f_matrix[i,j] =  prod(((exp.(fix_x_vector.*θ₀ .+ a_i_matrix[j,:] .+ a_j_matrix[j,:])).^total_y_grid_matrix[i,:])./( 1 .+ exp.(fix_x_vector.*θ₀ .+ a_i_matrix[j,:] .+ a_j_matrix[j,:]))) 
            end
        end

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

    return result_l_y, result_u_y

end

"""
Function that takes bounds for observed y's 
"""
function bounds_observedy(fix_y_matrix, total_y_grid, result_l_y, result_u_y)
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

    return result_l_y_real, result_u_y_real, result_l_y_mean, result_u_y_mean

end

"""
Function that computes true average marginal effects
"""
function true_average_marginal(Y, Aᵢ, Aⱼ, θ₀, N)
    Y_lower = tril(Y,-1)
    Aᵢ_lower = tril(Aᵢ,-1)
    Aⱼ_lower = tril(Aⱼ, -1)

    x₁_lower = ones(size(Y)[1], size(Y)[2])
    x₁_lower = tril(x₁_lower,-1)
    x₂_lower = -1 .*ones(size(Y)[1], size(Y)[2])
    x₂_lower = tril(x₂_lower,-1)

    individual_marginal = exp.(x₁_lower.*θ₀ .+ Aᵢ_lower .+ Aⱼ_lower)./( 1 .+ exp.(x₁_lower.*θ₀ .+ Aᵢ_lower .+ Aⱼ_lower)) .- exp.(x₂_lower.*θ₀ .+ Aᵢ_lower .+ Aⱼ_lower)./( 1 .+ exp.(x₂_lower.*θ₀ .+ Aᵢ_lower .+ Aⱼ_lower))
    true_marginal = sum(individual_marginal)/((N*(N-1))/2)

    return true_marginal

end
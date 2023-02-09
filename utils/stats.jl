using JuMP

"""
Function that creates a grids
"""
function agrids(αₗ, αₕ, λ₀, λ₁, grid_points, subnet_size)
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

    return a_i_matrix, a_j_matrix, total_a_grid

end

"""
Function that creates y grids
"""
function ygrids(subnet_size)
    y_possible = (0,1)
    total_y_grid = collect(Iterators.product(y_possible,y_possible,y_possible,y_possible,y_possible,y_possible))
    total_y_grid = reshape(total_y_grid,Int(2^(subnet_size*(subnet_size-1)/2)),1)
    total_y_grid_matrix = reshape(vec(reduce(hcat, getindex.(total_y_grid,i) for i in eachindex(total_y_grid[1]))),64,6)

    return total_y_grid, total_y_grid_matrix

end

"""
Function that creates combinations that satisfy Sijkl = 1 or Sijkl = -1
"""
function combination_condition(N, Y)
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
    
    combinations = Array{Set}(undef,length(combination_1))
    for i in 1:length(combination_1)
        combinations[i] = Set([combination_1[i],combination_2[i],combination_3[i],combination_4[i]])
    end

    return combinations

end

"""
Function that takes disjoint combinations
"""
function disjointcombinations(combinations)
    combinations_base = combinations
    combinations_final = Array{Set}(undef, 1)

    #initialize with one random set
    set_index = rand(1:length(combinations_base))
    combinations_final[1] = combinations_base[set_index]
    combinations_base = deleteat!(combinations_base,set_index)
    # then look if the remaining sets (in combinations base) are disjoint from the set in combinations_final[1]
    comparing_sets = [combinations_final[1] for _ in 1:length(combinations_base)]
    disjoint = isdisjoint.(combinations_base, comparing_sets)
    #if not disjoint, delete sets from combinations_base
    combinations_base = combinations_base[disjoint]

    condition = 0

    while condition == 0
        #1. pick a random set for the remaining combinations_base => we already know that it is disjoint with the previous set in combinations_final
        set_index = rand(1:length(combinations_base))
        new_set = combinations_base[set_index]
        #2. add to collection of combinations_final
        push!(combinations_final,new_set)
        #3. remove non disjoint sets from this one (the new chosen set) from combinations_base
        combinations_base = deleteat!(combinations_base,set_index)
        comparing_sets = [new_set for _ in 1:length(combinations_base)]
        disjoint = isdisjoint.(combinations_base, comparing_sets)
        combinations_base = combinations_base[disjoint]
        #4. check if combinations_base is not empty. if not, continue from the top.
        condition = isempty(combinations_base)
    end

    return combinations_final

end

"""
Function that creates matrices of Y and X for subnetworks defined by Graham condition
"""
function fix_yx_graham(combinations_final, Y, X)
    # Taking x vector as fixed : x_12, x_13, x_14, x_23, x_24, x_34
    fix_x_matrix = Array{NTuple{6,Float64}}(undef, length(combinations_final))
    fix_y_matrix = Array{NTuple{6,Float64}}(undef, length(combinations_final))

    for i in 1:length(combinations_final)
        ind1 = collect(combinations_final[i])[1]
        ind2 = collect(combinations_final[i])[2]
        ind3 = collect(combinations_final[i])[3]
        ind4 = collect(combinations_final[i])[4]
        fix_x_matrix[i] = (X[ind1, ind2], X[ind1, ind3], X[ind1, ind4], X[ind2, ind3], X[ind2, ind4], X[ind3, ind4])
        fix_y_matrix[i] = (Y[ind1, ind2], Y[ind1, ind3], Y[ind1, ind4], Y[ind2, ind3], Y[ind2, ind4], Y[ind3, ind4])
    end

    return fix_x_matrix, fix_y_matrix

end

"""
Function that creates matrices of Y and X for random subnetworks 
"""
function fix_yx(subnetworks, Y, X, s)
    # Taking x vector as fixed : x_12, x_13, x_14, x_23, x_24, x_34
    fix_x_matrix = Array{NTuple{6,Float64}}(undef, s)
    fix_y_matrix = Array{NTuple{6,Float64}}(undef, s)
    for i in 1:Int(s)
        ind1 = subnetworks[i][1]
        ind2 = subnetworks[i][2]
        ind3 = subnetworks[i][3]
        ind4 = subnetworks[i][4]
        fix_x_matrix[i] = (X[ind1, ind2], X[ind1, ind3], X[ind1, ind4], X[ind2, ind3], X[ind2, ind4], X[ind3, ind4])
        fix_y_matrix[i] = (Y[ind1, ind2], Y[ind1, ind3], Y[ind1, ind4], Y[ind2, ind3], Y[ind2, ind4], Y[ind3, ind4])
    end

    return fix_x_matrix, fix_y_matrix

end

"""
Function that computes the value of m for the first condition (given an x₁ and x₂)
"""
function mvalue(x₁_value, x₂_value, total_a_grid, θ₀, a_i_matrix, a_j_matrix)
    x₁ = x₁_value.*ones(length(total_a_grid), 6)
    x₂ = x₂_value.*ones(length(total_a_grid), 6)

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
Function that makes optimization
"""
function mainopt_graham(total_y_grid, combinations_final, fix_x_matrix, total_a_grid, θ₀, a_i_matrix, a_j_matrix, total_y_grid_matrix, bmax, bmin, m)

    #result_l_y = zeros(length(total_y_grid), length(combinations_final))
    result_l_y = Array{Union{Float64,Missing}}(missing,length(total_y_grid),length(combinations_final)) 
    #result_u_y = zeros(length(total_y_grid), length(combinations_final))
    result_u_y = Array{Union{Float64,Missing}}(missing,length(total_y_grid),length(combinations_final)) 

    # Fixing one subnetwork: one optimization for each subnetwork
    for sub in 1:length(combinations_final)
    
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
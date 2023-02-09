using Distributed

addprocs(7; exeflags="--project")

@everywhere begin
    using Random
    using Combinatorics
    using Distributions
    using Base.Iterators
    using JuMP
    using Clp
    using HiGHS
    using LinearAlgebra
    using DelimitedFiles
    using ProgressMeter
    using ThreadsX
    include("utils/dgp_graham.jl")
    include("utils/dgp_graham_undirected.jl")
    include("utils/stats_subnet_5.jl")
end

@everywhere function simulation(grid_points, N, θ₀, αₗ, αₕ, λ₀, λ₁, x₁_value, x₂_value, subnet_size, bmax, bmin)

    ### Take simulation
    Y, X, Aᵢ, Aⱼ = dgp_graham_undirected(N, θ₀, αₗ, αₕ, λ₀, λ₁)

    ### Grids of possible values of a and y
    # Create grids for aᵢ
    a_i_matrix, a_j_matrix, total_a_grid = agrids(αₗ, αₕ, λ₀, λ₁, grid_points, subnet_size)
    # as we are now working with undirected networks, take only the upper triangular matrix of the dgp
    # Create combinations (all possible) of y's : y_12, y_13, y_14, y_23, y_24, y_34, that is, subnetsize*(subnetsize-1)/2
    total_y_grid, total_y_grid_matrix = ygrids(subnet_size)

    # Subset subnetworks
    nodes = shuffle(1:N)
    s = Int(N/subnet_size)
    subnetworks = collect(partition(nodes,subnet_size))

    # Take matrices of Y's and X's that reflect the subnetworks of combinations: FUNCTION IS DIFFERENT IF WE DO NO CONSIDER GRAHAM'S TETRADS
    fix_x_matrix, fix_y_matrix = fix_yx(subnetworks, Y, X, s)

    ### Optimization
    # For first condition: get values of m
    m = mvalue(x₁_value, x₂_value, total_a_grid, θ₀, a_i_matrix, a_j_matrix)
    # For second condition: get the values of bmax, bmin in the beginning
    result_l_y, result_u_y = mainopt(total_y_grid, subnetworks, s, fix_x_matrix, total_a_grid, θ₀, a_i_matrix, a_j_matrix, total_y_grid_matrix, bmax, bmin, m)

    ### Take bounds for observed y's
    result_l_y_real, result_u_y_real, result_l_y_mean, result_u_y_mean = bounds_observedy(fix_y_matrix, total_y_grid, result_l_y, result_u_y)
    
    ### Compute true average marginal effects
    true_marginal = true_average_marginal(Y, Aᵢ, Aⱼ, θ₀, N)

    ### Extra results
    foreach(i -> Y[i, i] = 0, 1:N)
    number_links = sum(Y)/2
    percentage_links = number_links/((N*(N-1))/2) 
    number_links_subnetwork = sum.(fix_y_matrix)
    average_links_subnetwork = mean(sum.(fix_y_matrix))

    return result_l_y_mean, result_u_y_mean, true_marginal, s, percentage_links, average_links_subnetwork

end


grid_points = 5
N = 100
αₗ = -1
αₕ = -1
λ₀ = 1
λ₁ = 1
θ₀ = 1
x₁_value = 1
x₂_value = -1
subnet_size = 5
bmax = 1
bmin = -1

sims = 400
result = @time @showprogress pmap(1:sims) do sim
    simulation(grid_points, N, θ₀, αₗ, αₕ, λ₀, λ₁, x₁_value, x₂_value, subnet_size, bmax, bmin)
end

writedlm("results/outN$(N)al$(αₗ)ah$(αₕ)lambda0$(λ₀)lambda1$(λ₁)sims$(sims).csv", result, ',')




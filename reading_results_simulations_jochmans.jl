using FileIO, JLD2
using Statistics

## FIRST DESIGN JOCHMANS, N=30

results1 = FileIO.load("results_simulations/results_jochmans_DGP1_N30.jld2","result1_30")
results1 = reduce(hcat, getindex.(results1,i) for i in eachindex(results1[1]))
results1 = reduce(hcat,results1)
results1 = reshape(results1,1000,8)
replace!(results1, Inf=>NaN)

#now isnan, and isinf works

bias_β₁_comb1 = mean(results1[:,1])-1
bias_β₁_perm1 = mean(results1[:,2])-1
std_β₁_comb1 = std(results1[:,1])
std_β₁_perm1 = std(results1[:,2])

se₁1 = results1[:,3]
mean_se₁1 = mean(se₁1[isnan.(se₁1).==0])
perc_se₁1 = sum(isnan.(se₁1))/length(se₁1)

β₁ = results1[:,1]
β₁_se₁ = β₁[isnan.(se₁1).==0]
t_test₁1 = (β₁_se₁ .- 1) ./ (se₁1[isnan.(se₁1).==0])
size_t_test₁1 = sum(abs.(t_test₁1) .> 1.96)/length(t_test₁1)

se₂1 = results1[:,4]
mean_se₂1 = mean(se₂1[isnan.(se₂1).==0])
perc_se₂1 = sum(isnan.(se₂1))/length(se₂1)

β₁_se₂ = β₁[isnan.(se₂1).==0]
t_test₂1 = (β₁_se₂ .- 1) ./ (se₂1[isnan.(se₂1).==0])
size_t_test₂1 = sum(abs.(t_test₂1) .> 1.96)/length(t_test₂1)

se₃1 = results1[:,5]
mean_se₃1 = mean(se₃1[isnan.(se₃1).==0])
perc_se₃1 = sum(isnan.(se₃1))/length(se₃1)

β₁_se₃ = β₁[isnan.(se₃1).==0]
t_test₃1 = (β₁_se₃ .- 1) ./ (se₃1[isnan.(se₃1).==0])
size_t_test₃1 = sum(abs.(t_test₃1) .> 1.96)/length(t_test₃1)

se₄1 = results1[:,6]
mean_se₄1 = mean(se₄1[isnan.(se₄1).==0])
perc_se₄1 = sum(isnan.(se₄1))/length(se₄1)

β₁_se₄ = β₁[isnan.(se₄1).==0]
t_test₄1 = (β₁_se₄ .- 1) ./ (se₄1[isnan.(se₄1).==0])
size_t_test₄1 = sum(abs.(t_test₄1) .> 1.96)/length(t_test₄1)

mean_perc_quadruples1 = mean(results1[:,7])

mean_perc_links1 = mean(results1[:,8])
npₙ1 = mean_perc_quadruples1 * 30

## Second DESIGN JOCHMANS, N=30

results2 = FileIO.load("results_simulations/results_jochmans_DGP2_N30.jld2","result2_30")
results2 = reduce(hcat, getindex.(results2,i) for i in eachindex(results2[1]))
results2 = reduce(hcat,results2)
results2 = reshape(results2,1000,8)
replace!(results2, Inf=>NaN)

#now isnan, and isinf works

bias_β₁_comb2 = mean(results2[:,1])-1
bias_β₁_perm2 = mean(results2[:,2])-1
std_β₁_comb2 = std(results2[:,1])
std_β₁_perm2 = std(results2[:,2])

se₁2 = results2[:,3]
mean_se₁2 = mean(se₁2[isnan.(se₁2).==0])
perc_se₁2 = sum(isnan.(se₁2))/length(se₁2)

β₁ = results2[:,1]
β₁_se₁ = β₁[isnan.(se₁2).==0]
t_test₁2 = (β₁_se₁ .- 1) ./ (se₁2[isnan.(se₁2).==0])
size_t_test₁2 = sum(abs.(t_test₁2) .> 1.96)/length(t_test₁2)

se₂2 = results2[:,4]
mean_se₂2 = mean(se₂2[isnan.(se₂2).==0])
perc_se₂2 = sum(isnan.(se₂2))/length(se₂2)

β₁_se₂ = β₁[isnan.(se₂2).==0]
t_test₂2 = (β₁_se₂ .- 1) ./ (se₂2[isnan.(se₂2).==0])
size_t_test₂2 = sum(abs.(t_test₂2) .> 1.96)/length(t_test₂2)

se₃2 = results2[:,5]
mean_se₃2 = mean(se₃2[isnan.(se₃2).==0])
perc_se₃2 = sum(isnan.(se₃2))/length(se₃2)

β₁_se₃ = β₁[isnan.(se₃2).==0]
t_test₃2 = (β₁_se₃ .- 1) ./ (se₃2[isnan.(se₃2).==0])
size_t_test₃2 = sum(abs.(t_test₃2) .> 1.96)/length(t_test₃2)

se₄2 = results2[:,6]
mean_se₄2 = mean(se₄2[isnan.(se₄2).==0])
perc_se₄2 = sum(isnan.(se₄2))/length(se₄2)

β₁_se₄ = β₁[isnan.(se₄2).==0]
t_test₄2 = (β₁_se₄ .- 1) ./ (se₄2[isnan.(se₄2).==0])
size_t_test₄2 = sum(abs.(t_test₄2) .> 1.96)/length(t_test₄2)

mean_perc_quadruples2 = mean(results2[:,7])

mean_perc_links2 = mean(results2[:,8])
npₙ2 = mean_perc_quadruples2 * 30


## FIRST DESIGN JOCHMANS, N=30

results3 = FileIO.load("results_simulations/results_jochmans_DGP3_N30.jld2","result3_30")
results3 = reduce(hcat, getindex.(results3,i) for i in eachindex(results3[1]))
results3 = reduce(hcat,results3)
results3 = reshape(results3,1000,8)
replace!(results3, Inf=>NaN)

#now isnan, and isinf works

bias_β₁_comb3 = mean(results3[:,1])-1
bias_β₁_perm3 = mean(results3[:,2])-1
std_β₁_comb3 = std(results3[:,1])
std_β₁_perm3 = std(results3[:,2])

se₁3 = results3[:,3]
mean_se₁3 = mean(se₁3[isnan.(se₁3).==0])
perc_se₁3 = sum(isnan.(se₁3))/length(se₁3)

β₁ = results3[:,1]
β₁_se₁ = β₁[isnan.(se₁3).==0]
t_test₁3 = (β₁_se₁ .- 1) ./ (se₁3[isnan.(se₁3).==0])
size_t_test₁3 = sum(abs.(t_test₁3) .> 1.96)/length(t_test₁3)

se₂3 = results3[:,4]
mean_se₂3 = mean(se₂3[isnan.(se₂3).==0])
perc_se₂3 = sum(isnan.(se₂3))/length(se₂3)

β₁_se₂ = β₁[isnan.(se₂3).==0]
t_test₂3 = (β₁_se₂ .- 1) ./ (se₂3[isnan.(se₂3).==0])
size_t_test₂3 = sum(abs.(t_test₂3) .> 1.96)/length(t_test₂3)

se₃3 = results3[:,5]
mean_se₃3 = mean(se₃3[isnan.(se₃3).==0])
perc_se₃3 = sum(isnan.(se₃3))/length(se₃3)

β₁_se₃ = β₁[isnan.(se₃3).==0]
t_test₃3 = (β₁_se₃ .- 1) ./ (se₃3[isnan.(se₃3).==0])
size_t_test₃3 = sum(abs.(t_test₃3) .> 1.96)/length(t_test₃3)

se₄3 = results3[:,6]
mean_se₄3 = mean(se₄3[isnan.(se₄3).==0])
perc_se₄3 = sum(isnan.(se₄3))/length(se₄3)

β₁_se₄ = β₁[isnan.(se₄3).==0]
t_test₄3 = (β₁_se₄ .- 1) ./ (se₄3[isnan.(se₄3).==0])
size_t_test₄3 = sum(abs.(t_test₄3) .> 1.96)/length(t_test₄3)

mean_perc_quadruples3 = mean(results3[:,7])

mean_perc_links3 = mean(results3[:,8])
npₙ3 = mean_perc_quadruples3 * 30

## fourth DESIGN JOCHMANS, N=30

results4 = FileIO.load("results_simulations/results_jochmans_DGP4_N30.jld2","result4_30")
results4 = reduce(hcat, getindex.(results4,i) for i in eachindex(results4[1]))
results4 = reduce(hcat,results4)
results4 = reshape(results4,1000,8)
replace!(results4, Inf=>NaN)

#now isnan, and isinf works

bias_β₁_comb4 = mean(results4[:,1])-1
bias_β₁_perm4 = mean(results4[:,2])-1
std_β₁_comb4 = std(results4[:,1])
std_β₁_perm4 = std(results4[:,2])

se₁4 = results4[:,3]
mean_se₁4 = mean(se₁4[isnan.(se₁4).==0])
perc_se₁4 = sum(isnan.(se₁4))/length(se₁4)

β₁ = results4[:,1]
β₁_se₁ = β₁[isnan.(se₁4).==0]
t_test₁4 = (β₁_se₁ .- 1) ./ (se₁4[isnan.(se₁4).==0])
size_t_test₁4 = sum(abs.(t_test₁4) .> 1.96)/length(t_test₁4)

se₂4 = results4[:,4]
mean_se₂4 = mean(se₂4[isnan.(se₂4).==0])
perc_se₂4 = sum(isnan.(se₂4))/length(se₂4)

β₁_se₂ = β₁[isnan.(se₂4).==0]
t_test₂4 = (β₁_se₂ .- 1) ./ (se₂4[isnan.(se₂4).==0])
size_t_test₂4 = sum(abs.(t_test₂4) .> 1.96)/length(t_test₂4)

se₃4 = results4[:,5]
mean_se₃4 = mean(se₃4[isnan.(se₃4).==0])
perc_se₃4 = sum(isnan.(se₃4))/length(se₃4)

β₁_se₃ = β₁[isnan.(se₃4).==0]
t_test₃4 = (β₁_se₃ .- 1) ./ (se₃4[isnan.(se₃4).==0])
size_t_test₃4 = sum(abs.(t_test₃4) .> 1.96)/length(t_test₃4)

se₄4 = results4[:,6]
mean_se₄4 = mean(se₄4[isnan.(se₄4).==0])
perc_se₄4 = sum(isnan.(se₄4))/length(se₄4)

β₁_se₄ = β₁[isnan.(se₄4).==0]
t_test₄4 = (β₁_se₄ .- 1) ./ (se₄4[isnan.(se₄4).==0])
size_t_test₄4 = sum(abs.(t_test₄4) .> 1.96)/length(t_test₄4)

mean_perc_quadruples4 = mean(results4[:,7])

mean_perc_links4 = mean(results4[:,8])
npₙ4 = mean_perc_quadruples4 * 30

#### TABLE

Title = ["bias(β̂₁_comb)", "std(β̂₁_comb)", "bias(β̂₁_perm)", "std(β̂₁_perm)", "average(se₁(β̂₁))", "size(se₁(β̂₁))", "percentage(se₁(β̂₁))", "average(se₂(β̂₁))", "size(se₂(β̂₁))", "percentage(se₂(β̂₁))","average(se₃(β̂₁))", "size(se₃(β̂₁))", "percentage(se₃(β̂₁))","average(se₄(β̂₁))", "size(se₄(β̂₁))", "percentage(se₄(β̂₁))", "avg. percentage of links", "avg. percentage combinations", "avg. npₙ"]
Between = ["&&", "&&", "&&", "&&", "&&", "&&", "&&", "&&", "&&", "&&", "&&", "&&", "&&", "&&", "&&", "&&", "&&", "&&", "&&"]
End = ["\\","\\","\\","\\","\\","\\","\\","\\","\\","\\", "\\","\\","\\","\\","\\","\\","\\","\\","\\"]
Design1 = [bias_β₁_comb1, std_β₁_comb1, bias_β₁_perm1, std_β₁_perm1, mean_se₁1, size_t_test₁1, perc_se₁1, mean_se₂1, size_t_test₂1, perc_se₂1, mean_se₃1, size_t_test₃1, perc_se₃1, mean_se₄1, size_t_test₄1, perc_se₄1, mean_perc_links1, mean_perc_quadruples1, npₙ1] 
Design2 = [bias_β₁_comb2, std_β₁_comb2, bias_β₁_perm2, std_β₁_perm2, mean_se₁2, size_t_test₁2, perc_se₁2, mean_se₂2, size_t_test₂2, perc_se₂2, mean_se₃2, size_t_test₃2, perc_se₃2, mean_se₄2, size_t_test₄2, perc_se₄2, mean_perc_links2, mean_perc_quadruples2, npₙ2] 
Design3 = [bias_β₁_comb3, std_β₁_comb3, bias_β₁_perm3, std_β₁_perm3, mean_se₁3, size_t_test₁3, perc_se₁3, mean_se₂3, size_t_test₂3, perc_se₂3, mean_se₃3, size_t_test₃3, perc_se₃3, mean_se₄3, size_t_test₄3, perc_se₄3, mean_perc_links3, mean_perc_quadruples3, npₙ3] 
Design4 = [bias_β₁_comb4, std_β₁_comb4, bias_β₁_perm4, std_β₁_perm4, mean_se₁4, size_t_test₁4, perc_se₁4, mean_se₂4, size_t_test₂4, perc_se₂4, mean_se₃4, size_t_test₃4, perc_se₃4, mean_se₄4, size_t_test₄4, perc_se₄4, mean_perc_links4, mean_perc_quadruples4, npₙ4] 

Table_results_Jochmans_N30 = [Title Between Design1 Between Design2 Between Design3 Between Design4 End]





## FIRST DESIGN JOCHMANS, N=45

results1 = FileIO.load("results_simulations/results_jochmans_DGP1_N45.jld2","result1_45")
results1 = reduce(hcat, getindex.(results1,i) for i in eachindex(results1[1]))
results1 = reduce(hcat,results1)
results1 = reshape(results1,1000,8)
replace!(results1, Inf=>NaN)

#now isnan, and isinf works

bias_β₁_comb1 = mean(results1[:,1])-1
bias_β₁_perm1 = mean(results1[:,2])-1
std_β₁_comb1 = std(results1[:,1])
std_β₁_perm1 = std(results1[:,2])

se₁1 = results1[:,3]
mean_se₁1 = mean(se₁1[isnan.(se₁1).==0])
perc_se₁1 = sum(isnan.(se₁1))/length(se₁1)

β₁ = results1[:,1]
β₁_se₁ = β₁[isnan.(se₁1).==0]
t_test₁1 = (β₁_se₁ .- 1) ./ (se₁1[isnan.(se₁1).==0])
size_t_test₁1 = sum(abs.(t_test₁1) .> 1.96)/length(t_test₁1)

se₂1 = results1[:,4]
mean_se₂1 = mean(se₂1[isnan.(se₂1).==0])
perc_se₂1 = sum(isnan.(se₂1))/length(se₂1)

β₁_se₂ = β₁[isnan.(se₂1).==0]
t_test₂1 = (β₁_se₂ .- 1) ./ (se₂1[isnan.(se₂1).==0])
size_t_test₂1 = sum(abs.(t_test₂1) .> 1.96)/length(t_test₂1)

se₃1 = results1[:,5]
mean_se₃1 = mean(se₃1[isnan.(se₃1).==0])
perc_se₃1 = sum(isnan.(se₃1))/length(se₃1)

β₁_se₃ = β₁[isnan.(se₃1).==0]
t_test₃1 = (β₁_se₃ .- 1) ./ (se₃1[isnan.(se₃1).==0])
size_t_test₃1 = sum(abs.(t_test₃1) .> 1.96)/length(t_test₃1)

se₄1 = results1[:,6]
mean_se₄1 = mean(se₄1[isnan.(se₄1).==0])
perc_se₄1 = sum(isnan.(se₄1))/length(se₄1)

β₁_se₄ = β₁[isnan.(se₄1).==0]
t_test₄1 = (β₁_se₄ .- 1) ./ (se₄1[isnan.(se₄1).==0])
size_t_test₄1 = sum(abs.(t_test₄1) .> 1.96)/length(t_test₄1)

mean_perc_quadruples1 = mean(results1[:,7])

mean_perc_links1 = mean(results1[:,8])
npₙ1 = mean_perc_quadruples1 * 45

## Second DESIGN JOCHMANS, N=45

results2 = FileIO.load("results_simulations/results_jochmans_DGP2_N45.jld2","result2_45")
results2 = reduce(hcat, getindex.(results2,i) for i in eachindex(results2[1]))
results2 = reduce(hcat,results2)
results2 = reshape(results2,1000,8)
replace!(results2, Inf=>NaN)

#now isnan, and isinf works

bias_β₁_comb2 = mean(results2[:,1])-1
bias_β₁_perm2 = mean(results2[:,2])-1
std_β₁_comb2 = std(results2[:,1])
std_β₁_perm2 = std(results2[:,2])

se₁2 = results2[:,3]
mean_se₁2 = mean(se₁2[isnan.(se₁2).==0])
perc_se₁2 = sum(isnan.(se₁2))/length(se₁2)

β₁ = results2[:,1]
β₁_se₁ = β₁[isnan.(se₁2).==0]
t_test₁2 = (β₁_se₁ .- 1) ./ (se₁2[isnan.(se₁2).==0])
size_t_test₁2 = sum(abs.(t_test₁2) .> 1.96)/length(t_test₁2)

se₂2 = results2[:,4]
mean_se₂2 = mean(se₂2[isnan.(se₂2).==0])
perc_se₂2 = sum(isnan.(se₂2))/length(se₂2)

β₁_se₂ = β₁[isnan.(se₂2).==0]
t_test₂2 = (β₁_se₂ .- 1) ./ (se₂2[isnan.(se₂2).==0])
size_t_test₂2 = sum(abs.(t_test₂2) .> 1.96)/length(t_test₂2)

se₃2 = results2[:,5]
mean_se₃2 = mean(se₃2[isnan.(se₃2).==0])
perc_se₃2 = sum(isnan.(se₃2))/length(se₃2)

β₁_se₃ = β₁[isnan.(se₃2).==0]
t_test₃2 = (β₁_se₃ .- 1) ./ (se₃2[isnan.(se₃2).==0])
size_t_test₃2 = sum(abs.(t_test₃2) .> 1.96)/length(t_test₃2)

se₄2 = results2[:,6]
mean_se₄2 = mean(se₄2[isnan.(se₄2).==0])
perc_se₄2 = sum(isnan.(se₄2))/length(se₄2)

β₁_se₄ = β₁[isnan.(se₄2).==0]
t_test₄2 = (β₁_se₄ .- 1) ./ (se₄2[isnan.(se₄2).==0])
size_t_test₄2 = sum(abs.(t_test₄2) .> 1.96)/length(t_test₄2)

mean_perc_quadruples2 = mean(results2[:,7])

mean_perc_links2 = mean(results2[:,8])
npₙ2 = mean_perc_quadruples2 * 45


## FIRST DESIGN JOCHMANS, N=45

results3 = FileIO.load("results_simulations/results_jochmans_DGP3_N45.jld2","result3_45")
results3 = reduce(hcat, getindex.(results3,i) for i in eachindex(results3[1]))
results3 = reduce(hcat,results3)
results3 = reshape(results3,1000,8)
replace!(results3, Inf=>NaN)

#now isnan, and isinf works

bias_β₁_comb3 = mean(results3[:,1])-1
bias_β₁_perm3 = mean(results3[:,2])-1
std_β₁_comb3 = std(results3[:,1])
std_β₁_perm3 = std(results3[:,2])

se₁3 = results3[:,3]
mean_se₁3 = mean(se₁3[isnan.(se₁3).==0])
perc_se₁3 = sum(isnan.(se₁3))/length(se₁3)

β₁ = results3[:,1]
β₁_se₁ = β₁[isnan.(se₁3).==0]
t_test₁3 = (β₁_se₁ .- 1) ./ (se₁3[isnan.(se₁3).==0])
size_t_test₁3 = sum(abs.(t_test₁3) .> 1.96)/length(t_test₁3)

se₂3 = results3[:,4]
mean_se₂3 = mean(se₂3[isnan.(se₂3).==0])
perc_se₂3 = sum(isnan.(se₂3))/length(se₂3)

β₁_se₂ = β₁[isnan.(se₂3).==0]
t_test₂3 = (β₁_se₂ .- 1) ./ (se₂3[isnan.(se₂3).==0])
size_t_test₂3 = sum(abs.(t_test₂3) .> 1.96)/length(t_test₂3)

se₃3 = results3[:,5]
mean_se₃3 = mean(se₃3[isnan.(se₃3).==0])
perc_se₃3 = sum(isnan.(se₃3))/length(se₃3)

β₁_se₃ = β₁[isnan.(se₃3).==0]
t_test₃3 = (β₁_se₃ .- 1) ./ (se₃3[isnan.(se₃3).==0])
size_t_test₃3 = sum(abs.(t_test₃3) .> 1.96)/length(t_test₃3)

se₄3 = results3[:,6]
mean_se₄3 = mean(se₄3[isnan.(se₄3).==0])
perc_se₄3 = sum(isnan.(se₄3))/length(se₄3)

β₁_se₄ = β₁[isnan.(se₄3).==0]
t_test₄3 = (β₁_se₄ .- 1) ./ (se₄3[isnan.(se₄3).==0])
size_t_test₄3 = sum(abs.(t_test₄3) .> 1.96)/length(t_test₄3)

mean_perc_quadruples3 = mean(results3[:,7])

mean_perc_links3 = mean(results3[:,8])
npₙ3 = mean_perc_quadruples3 * 45

## fourth DESIGN JOCHMANS, N=45

results4 = FileIO.load("results_simulations/results_jochmans_DGP4_N45.jld2","result4_45")
results4 = reduce(hcat, getindex.(results4,i) for i in eachindex(results4[1]))
results4 = reduce(hcat,results4)
results4 = reshape(results4,1000,8)
replace!(results4, Inf=>NaN)

#now isnan, and isinf works

bias_β₁_comb4 = mean(results4[:,1])-1
bias_β₁_perm4 = mean(results4[:,2])-1
std_β₁_comb4 = std(results4[:,1])
std_β₁_perm4 = std(results4[:,2])

se₁4 = results4[:,3]
mean_se₁4 = mean(se₁4[isnan.(se₁4).==0])
perc_se₁4 = sum(isnan.(se₁4))/length(se₁4)

β₁ = results4[:,1]
β₁_se₁ = β₁[isnan.(se₁4).==0]
t_test₁4 = (β₁_se₁ .- 1) ./ (se₁4[isnan.(se₁4).==0])
size_t_test₁4 = sum(abs.(t_test₁4) .> 1.96)/length(t_test₁4)

se₂4 = results4[:,4]
mean_se₂4 = mean(se₂4[isnan.(se₂4).==0])
perc_se₂4 = sum(isnan.(se₂4))/length(se₂4)

β₁_se₂ = β₁[isnan.(se₂4).==0]
t_test₂4 = (β₁_se₂ .- 1) ./ (se₂4[isnan.(se₂4).==0])
size_t_test₂4 = sum(abs.(t_test₂4) .> 1.96)/length(t_test₂4)

se₃4 = results4[:,5]
mean_se₃4 = mean(se₃4[isnan.(se₃4).==0])
perc_se₃4 = sum(isnan.(se₃4))/length(se₃4)

β₁_se₃ = β₁[isnan.(se₃4).==0]
t_test₃4 = (β₁_se₃ .- 1) ./ (se₃4[isnan.(se₃4).==0])
size_t_test₃4 = sum(abs.(t_test₃4) .> 1.96)/length(t_test₃4)

se₄4 = results4[:,6]
mean_se₄4 = mean(se₄4[isnan.(se₄4).==0])
perc_se₄4 = sum(isnan.(se₄4))/length(se₄4)

β₁_se₄ = β₁[isnan.(se₄4).==0]
t_test₄4 = (β₁_se₄ .- 1) ./ (se₄4[isnan.(se₄4).==0])
size_t_test₄4 = sum(abs.(t_test₄4) .> 1.96)/length(t_test₄4)

mean_perc_quadruples4 = mean(results4[:,7])

mean_perc_links4 = mean(results4[:,8])
npₙ4 = mean_perc_quadruples4 * 45

#### TABLE

Title = ["bias(β̂₁_comb)", "std(β̂₁_comb)", "bias(β̂₁_perm)", "std(β̂₁_perm)", "average(se₁(β̂₁))", "size(se₁(β̂₁))", "percentage(se₁(β̂₁))", "average(se₂(β̂₁))", "size(se₂(β̂₁))", "percentage(se₂(β̂₁))","average(se₃(β̂₁))", "size(se₃(β̂₁))", "percentage(se₃(β̂₁))","average(se₄(β̂₁))", "size(se₄(β̂₁))", "percentage(se₄(β̂₁))", "avg. percentage of links", "avg. percentage combinations", "avg. npₙ"]
Between = ["&&", "&&", "&&", "&&", "&&", "&&", "&&", "&&", "&&", "&&", "&&", "&&", "&&", "&&", "&&", "&&", "&&", "&&", "&&"]
End = ["\\","\\","\\","\\","\\","\\","\\","\\","\\","\\", "\\","\\","\\","\\","\\","\\","\\","\\","\\"]
Design1 = [bias_β₁_comb1, std_β₁_comb1, bias_β₁_perm1, std_β₁_perm1, mean_se₁1, size_t_test₁1, perc_se₁1, mean_se₂1, size_t_test₂1, perc_se₂1, mean_se₃1, size_t_test₃1, perc_se₃1, mean_se₄1, size_t_test₄1, perc_se₄1, mean_perc_links1, mean_perc_quadruples1, npₙ1] 
Design2 = [bias_β₁_comb2, std_β₁_comb2, bias_β₁_perm2, std_β₁_perm2, mean_se₁2, size_t_test₁2, perc_se₁2, mean_se₂2, size_t_test₂2, perc_se₂2, mean_se₃2, size_t_test₃2, perc_se₃2, mean_se₄2, size_t_test₄2, perc_se₄2, mean_perc_links2, mean_perc_quadruples2, npₙ2] 
Design3 = [bias_β₁_comb3, std_β₁_comb3, bias_β₁_perm3, std_β₁_perm3, mean_se₁3, size_t_test₁3, perc_se₁3, mean_se₂3, size_t_test₂3, perc_se₂3, mean_se₃3, size_t_test₃3, perc_se₃3, mean_se₄3, size_t_test₄3, perc_se₄3, mean_perc_links3, mean_perc_quadruples3, npₙ3] 
Design4 = [bias_β₁_comb4, std_β₁_comb4, bias_β₁_perm4, std_β₁_perm4, mean_se₁4, size_t_test₁4, perc_se₁4, mean_se₂4, size_t_test₂4, perc_se₂4, mean_se₃4, size_t_test₃4, perc_se₃4, mean_se₄4, size_t_test₄4, perc_se₄4, mean_perc_links4, mean_perc_quadruples4, npₙ4] 

Table_results_Jochmans_N45 = [Title Between Design1 Between Design2 Between Design3 Between Design4 End]


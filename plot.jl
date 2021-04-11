using Printf
using Distributions
using Plots
using Statistics
using DelimitedFiles

function myqqplot(obs,F⁰,title)
    nobs=length(obs)
    sort!(obs)
    quantiles⁰ = [quantile(F⁰,i/nobs) for i in 1:nobs]
    # Note that only n-1 points may be plotted, as quantile(F⁰,1) may be inf
    plot(quantiles⁰, obs, seriestype=:scatter, xlabel="Theoretical Quantiles", ylabel = "Sample Quantiles", title=title, label="" )
    plot!(obs,obs,label="")
end

#Looping over all simulation results to plot and generating tables

β₁=1

# for each design create a matrix with rows Nxsims and columns with main results
# for factors: S, N, meanF₁, F̄₁, meanF₂, F̄₂
factorsdesign = [zeros(3*4,6), zeros(3*4,6), zeros(3*4,6), zeros(3*4,6)]
# for betas: S, N, meanβ̂₁ - β₁, varβ̂₁, meanvarβ̂₁eff144, meanvarβ̂₁eff72, meanvarβ̂₁ineff144, meanvarβ̂₁ineff72
betasdesign = [zeros(3*4,6), zeros(3*4,6), zeros(3*4,6), zeros(3*4,6)]
# for sizes of t-tests: S, N, sizeβ̂₁eff144, sizeβ̂₁eff72, sizeβ̂₁ineff144, sizeβ̂₁ineff72
sizedesign = [zeros(3*4,4), zeros(3*4,4), zeros(3*4,4), zeros(3*4,4)]

#for design ∈ [1,2,3,4]
for design ∈ [1]
    counter = 1
    for sims ∈ [1000, 5000, 10000]
        for N ∈ [10, 20, 30, 50]

            #reading results
            result = readdlm("results/outN$(N)sims$(sims)_design$(design).csv", ',', Float64)

            #some further results to obtain the factors and plots
            δ₂ = result[:,1]
            Δ₂ = result[:,2]
            Ustat = result[:,3]
            β̂₁ = result[:,4]
            varβ̂₁144 = result[:,5]
            varβ̂₁72 = result[:,6]

            """
            Results for factor analysis
            """

            varUstat = var(Ustat)
            F₁ = varUstat * (N * (N - 1)) ./ δ₂
            F₂ = varUstat * (N * (N - 1)) ./ Δ₂

            #storing in matrices
            factorsdesign[design][counter,1] = sims
            factorsdesign[design][counter,2] = N
            factorsdesign[design][counter,3] = mean(F₁) #mean(F₁)
            factorsdesign[design][counter,4] = varUstat * (N * (N - 1)) / mean(δ₂) #F̄₁
            factorsdesign[design][counter,5] = mean(F₂) #mean(F₂)
            factorsdesign[design][counter,6] = varUstat * (N * (N - 1)) / mean(Δ₂) #F̄₂

            #ploting histograms
            hist₁ = histogram(F₁, bins=:scott, title="For N = $N, and S = $sims", label="", xlabel="F₁",
            ylabel="Frequency")
            savefig("results/histograms-F1/histF1_N$(N)_S$(sims)_design$(design).png")
            hist₂ = histogram(F₂, bins=:scott, title="For N = $N, and S = $sims", label="", xlabel="F₂",
            ylabel="Frequency")
            savefig("results/histograms-F2/histF2_N$(N)_S$(sims)_design$(design).png")

            """
            Results for OLS estimator
            """

            #storing in matrices
            betasdesign[design][counter,1] = sims
            betasdesign[design][counter,2] = N
            betasdesign[design][counter,3] = mean(β̂₁) - β₁ #biasβ̂₁
            betasdesign[design][counter,4] = var(β̂₁)
            betasdesign[design][counter,5] = mean(varβ̂₁144)
            betasdesign[design][counter,6] = mean(varβ̂₁72)

            #ploting histogram and qq-plot
            hist₃ = histogram(β̂₁, bins=:scott, title="For N = $N, and S = $sims", label="", xlabel="β̂₁",
            ylabel="Frequency")
            savefig("results/histograms-beta/histbeta_N$(N)_S$(sims)_design$(design).png")

            myqqplot(β̂₁,Normal(mean(β̂₁), std(β̂₁)),"For N = $N, and S = $sims")
            savefig("results/qqplot-beta/qqplotbeta_N$(N)_S$(sims)_design$(design).png")

            """
            Results for size of t-tests
            """

            #storing in matrices
            sizedesign[design][counter,1] = sims
            sizedesign[design][counter,2] = N
            sizedesign[design][counter,3] = sum(abs.((β̂₁ .- β₁)./(sqrt.(varβ̂₁144))) .> 1.96)/sims #sizeβ̂₁eff144
            sizedesign[design][counter,4] = sum(abs.((β̂₁ .- β₁)./(sqrt.(varβ̂₁72))) .> 1.96)/sims #sizeβ̂₁eff72

            #plotting QQ-plots of t-tests
            ttest144 = (β̂₁ .- β₁)./(sqrt.(varβ̂₁144))
            myqqplot(ttest144,Normal(mean(ttest144), std(ttest144)),"For N = $N, and S = $sims")
            savefig("results/qqplot-ttest/qqplotttest144_N$(N)_S$(sims)_design$(design).png")
            
            ttest72 = (β̂₁ .- β₁)./(sqrt.(varβ̂₁72))
            myqqplot(ttest72,Normal(mean(ttest72), std(ttest72)),"For N = $N, and S = $sims")
            savefig("results/qqplot-ttest/qqplotttest72_N$(N)_S$(sims)_design$(design).png")

            counter += 1
        end
    end
end
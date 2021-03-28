using Printf
using Distributions
using Plots
using Statistics
using DelimitedFiles

#Looping over all simulation results to plot and generating tables

β₁=1

# for each design create a matrix with rows Nxsims and columns with main results
# for factors: S, N, meanF₁, F̄₁, meanF₂, F̄₂
factorsdesign = [zeros(3*4+1,6), zeros(3*4+1,6), zeros(3*4+1,6), zeros(3*4+1,6)]
# for betas: S, N, meanβ̂₁ - β₁, varβ̂₁, meanvarβ̂₁eff144, meanvarβ̂₁eff72, meanvarβ̂₁ineff144, meanvarβ̂₁ineff72
betasdesign = [zeros(3*4+1,8), zeros(3*4+1,8), zeros(3*4+1,8), zeros(3*4+1,8)]
# for sizes of t-tests: S, N, sizeβ̂₁eff144, sizeβ̂₁eff72, sizeβ̂₁ineff144, sizeβ̂₁ineff72
sizedesign = [zeros(3*4+1,6), zeros(3*4+1,6), zeros(3*4+1,6), zeros(3*4+1,6)]

for design ∈ [1,2,3,4]
    counter = 1
    for sims ∈ [1000, 5000, 10000]
        for N ∈ [10, 20, 30, 50]

            #reading results
            result = readdlm("results/outN$(N)sims$(sims)_design$(design).csv", ',', Float64)

            #some further results to obtain the factors and plots
            Δ₂₁ = result[:,1]
            Δ₂₂ = result[:,2]
            Ustat = result[:,3]
            β̂₁ = result[:,4]
            varβ̂₁eff144 = result[:,5]
            varβ̂₁eff72 = result[:,6]
            varβ̂₁ineff144 = result[:,7]
            varβ̂₁ineff72 = result[:,8]

            """
            Results for factor analysis
            """

            varUstat = var(Ustat)
            F₁ = varUstat * (N * (N - 1)) ./ Δ₂₁
            meanF₁ = mean(F₁)
            F̄₁ = varUstat * (N * (N - 1)) / mean(Δ₂₁) # look at this one
            F₂ = varUstat * (N * (N - 1)) ./ Δ₂₂
            meanF₂ = mean(F₂)
            F̄₂ = varUstat * (N * (N - 1)) / mean(Δ₂₂)

            #storing in matrices
            factorsdesign[design][counter,1] = sims
            factorsdesign[design][counter,2] = N
            factorsdesign[design][counter,3] = F₁
            factorsdesign[design][counter,4] = F̄₁
            factorsdesign[design][counter,5] = F₂
            factorsdesign[design][counter,6] = F̄₂

            #ploting histograms
            hist₁ = histogram(F₁, bins=:scott, title="For N = $N, and S = $sims", label="", xlabel="F₁",
            ylabel="Frequency")
            savefig("results/histograms-F1/histF1_N$(N)_S$(sims)_design$(design).pdf")
            hist₂ = histogram(F₂, bins=:scott, title="For N = $N, and S = $sims", label="", xlabel="F₂",
            ylabel="Frequency")
            savefig("results/histograms-F2/histF2_N$(N)_S$(sims)_design$(design).pdf")

            """
            Results for OLS estimator
            """

            biasβ̂₁ = mean(β̂₁) - β₁
            varβ̂₁ = var(β̂₁)
            meanvarβ̂₁eff144 = mean(varβ̂₁eff144)
            meanvarβ̂₁eff72 = mean(varβ̂₁eff72)
            meanvarβ̂₁ineff144 = mean(varβ̂₁ineff144)
            meanvarβ̂₁ineff72 = mean(varβ̂₁ineff72)

            #storing in matrices
            betasdesign[design][counter,1] = sims
            betasdesign[design][counter,2] = N
            betasdesign[design][counter,3] = biasβ̂₁
            betasdesign[design][counter,4] = varβ̂₁
            betasdesign[design][counter,5] = meanvarβ̂₁eff144
            betasdesign[design][counter,6] = meanvarβ̂₁eff72
            betasdesign[design][counter,7] = meanvarβ̂₁ineff144
            betasdesign[design][counter,8] = meanvarβ̂₁ineff72

            #ploting histogram and qq-plot
            hist₃ = histogram(β̂₁, bins=:scott, title="For N = $N, and S = $sims", label="", xlabel="β̂₁",
            ylabel="Frequency")
            savefig("results/histograms-beta/histbeta_N$(N)_S$(sims)_design$(design).pdf")

            function myqqplot(obs,F⁰,title)
                nobs=length(obs)
                sort!(obs)
                quantiles⁰ = [quantile(F⁰,i/nobs) for i in 1:nobs]
                # Note that only n-1 points may be plotted, as quantile(F⁰,1) may be inf
                plot(quantiles⁰, obs, seriestype=:scatter, xlabel="Theoretical Quantiles", ylabel = "Sample Quantiles", title=title, label="" )
                plot!(obs,obs,label="")
            end
            myqqplot(β̂₁,Normal(mean(β̂₁), std(β̂₁)),"For N = $N, and S = $sims")
            savefig("results/qqplot-beta/qqplotbeta_N$(N)_S$(sims)_design$(design).pdf")


            """
            Results for size of t-tests
            """
         
            sizeβ̂₁eff144 = sum(abs.((β̂₁ .- β₁)./(sqrt.(varβ̂₁eff144))) .> 1.96)/sims 
            sizeβ̂₁eff72 = sum(abs.((β̂₁ .- β₁)./(sqrt.(varβ̂₁eff72))) .> 1.96)/sims 
            sizeβ̂₁ineff144 = sum(abs.((β̂₁ .- β₁)./(sqrt.(varβ̂₁ineff144))) .> 1.96)/sims 
            sizeβ̂₁ineff72  = sum(abs.((β̂₁ .- β₁)./(sqrt.(varβ̂₁ineff72))) .> 1.96)/sims 

            #storing in matrices
            sizedesign[design][counter,1] = sims
            sizedesign[design][counter,2] = N
            sizedesign[design][counter,3] = sizeβ̂₁eff144
            sizedesign[design][counter,4] = sizeβ̂₁eff72
            sizedesign[design][counter,5] = sizeβ̂₁ineff144
            sizedesign[design][counter,6] = Fsizeβ̂₁ineff72

            counter += 1
        end
    end
end
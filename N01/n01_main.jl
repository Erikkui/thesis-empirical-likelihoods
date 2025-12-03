using CairoMakie
using LoopVectorization
using Statistics
using LinearAlgebra
using ProgressBars
using Distributions
using Distances

include("wrapper.jl")
include("include_lib.jl")
INCLUDE_PATH = joinpath( pwd(), "lib" )
include_lib( INCLUDE_PATH )

colors = Makie.to_colormap(:seaborn_dark)
set_theme!(palette = (; color = colors))

data = Dict{Symbol, Any}()
Ndata = 100
nepo = 1
nbin = 10

nresample = 1
nsim = 100

LL = []
data[:LL] = copy(LL)

mu = 3.0
sigma = 1.0
params = [mu, sigma]
npar = length( params )

R_all = Matrix{Float64}(undef, nepo, Ndata)
mins = zeros( nepo )
maxs = zeros( nepo )
for ii in 1:nepo
    R = randn( 1, Ndata )*sigma .+ mu
    R_all[ii, :] = R
    mins[ii] = minimum( R )
    maxs[ii] = maximum( R )
end
m = maximum(mins)
M = minimum(maxs)
data[:R_all] = R_all

if !isempty(data[:LL])
    print("Calculating minmax for CIL_ID features...")
    RR = reshape( R_all[1, :], 1, Ndata )
    mins2 = zeros( nresample, length(data[:LL]) )
    maxs2 = zeros( nresample, length(data[:LL]) )
    LL = copy( data[:LL] )
    for rr in 1:nresample
        x_inds = rand( 1:Ndata, Ndata ) |> unique
        y_inds = setdiff( 1:Ndata, x_inds )
        x = RR[ 1, x_inds ]
        y = RR[ 1, y_inds ]
        D = pairwise(Euclidean(), x', y'; dims=2)

        cil_id = CIL_ID( D, LL )
        for feature in eachindex( cil_id )
            mins2[rr, feature] = minimum( cil_id[feature] )
            maxs2[rr, feature] = maximum( cil_id[feature] )
        end
    end
    mins2 = maximum( mins2, dims=1 )
    maxs2 = minimum( maxs2, dims=1 )
end


data[:minmax] = [m, M]
data[:r] = R_all
data[:nbin] = nbin
data[:uni] = "xax"
data[:theta] = params .* rand( length(params) )
data[:nsim] = nsim
data[:nresample] = nresample
data[:Ndata] = Ndata

bins = bin_select( data )
data[:bins] = [bins]
if !isempty(data[:LL])
    println("Calculating bins for CIL_ID features...")
    for ii in 1:size(mins2, 2)
        data[:minmax] = [mins2[1, ii], maxs2[1, ii]]
        bins = bin_select( data )
        push!( data[:bins], bins )
    end
end


# Create ecdfs from observations
# fig = Figure( size = (800, 400), padding = 1, fontsize = 13 )
# ax = Axis( fig[1, 1])
S_y = zeros( nepo, nbin )
S_y2 = zeros( nepo, length(data[:LL])*nbin )
for (ii, R) in enumerate( eachrow(R_all) )
    LL = copy( data[:LL] )
    cdf_y, _ = empcdf( R, nx=nbin, x=data[:bins][1] )
    S_y[ii, :] = vec( cdf_y )
    # lines!( ax, data[:bins][1], vec(cdf_y), color = colors[1] )
    if !isempty(data[:LL])
        println("Calculating CIL_ID ecdfs")
        D = pairwise(Euclidean(), reshape(R, :, 1), reshape(R, :, 1); dims=2)
        cil_id = CIL_ID( D, LL )
        cdf_f_ii = Vector{Float64}()
        for feature in eachindex( cil_id )
            cdf_feat, _ = empcdf( cil_id[feature], nx=nbin, x=data[:bins][1+feature] )
            cdf_f_ii = vcat( cdf_f_ii, vec( cdf_feat ) )
        end
        S_y2[ii, :] = cdf_f_ii
    end
end
if !isempty(data[:LL])
    S_y = hcat( S_y, S_y2 )
end
data[:muu_data] = mean( S_y, dims=1 )
# display( fig )


# MCMC settings
mcmc_options = Dict{Symbol, Any}(
        :nsimu => 20000,
        :update_int => 30,
        :adapt_int => 20,
        :qcov => Matrix{Float64}( I, npar, npar )*1e-2,
    )
mcmc_model = Dict{Symbol, Any}(
    :sigma2 => 1,
    :noise => 0.0,
    :ssfun => like_eval
)

chain, _, _ = mcmcrun2( mcmc_model, data, mcmc_options )
chain = chain[5001:end, :]  # discard burn-in

# Plot MCMC parameter chains
fig_chain = Figure(size = (800, 400), padding = 1, fontsize = 13)
for i in 1:npar
    ax = Axis(fig_chain[i, 1], title = "Parameter $(i) Chain")
    scatter!(ax, chain[:, i], color = colors[1])
end
display(fig_chain)

# Plot joint posterior
fig_posterior = Figure(size = (800, 400), padding = 1, fontsize = 13)
ax_posterior = Axis(fig_posterior[1, 1], title = "Joint Posterior", xlabel = "Mu", ylabel = "Sigma")
scatter!(ax_posterior, chain[:, 1], chain[:, 2], color = colors[1], markersize = 9)
display(fig_posterior)

# Calculate posterior mean and plot estimated vs true normal distributions
posterior_mean = mean(chain, dims = 1)
mu_est, sigma_est = posterior_mean[1], posterior_mean[2]

fig_dist = Figure(size = (800, 400), padding = 1, fontsize = 13)
ax_dist = Axis(fig_dist[1, 1], title = "Estimated vs True Normal Distribution")
x_vals = range(m - 0.1, M + 0.1, length = 1000)
true_pdf = pdf.(Normal(mu, sigma), x_vals)
est_pdf = pdf.(Normal(mu_est, sigma_est), x_vals)
lines!(ax_dist, x_vals, true_pdf, color = colors[1], label = "True")
lines!(ax_dist, x_vals, est_pdf, color = colors[2], label = "Estimated")
axislegend(ax_dist)
display(fig_dist)

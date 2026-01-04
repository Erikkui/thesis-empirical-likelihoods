function OU_solve( theta, x0, Ndata; dt = 0.1, rng_seed = nothing )

    if isnothing( rng_seed )
        local_rng = Random.GLOBAL_RNG
    else
        local_rng = MersenneTwister( rng_seed )
    end

    lambda, sigma = theta
    ndim = length(x0)

    tt = ( 0:Ndata ) .* dt
    noise = randn( local_rng, ndim, length(tt) )

    x = zeros( ndim, length(tt) )
    x[:, 1] .= x0

    for ii in 2:length(tt)
        eps_t = @view noise[:, ii]
        x_prev = @view x[:, ii-1]
        x[:, ii] = x_prev .- lambda*x_prev*dt .+ sigma*sqrt(dt)*eps_t
    end

    return x
end


## Testing OU
# using CairoMakie
# include("lib/calculate_diffs.jl")
# theta = [5, 8]
# x0 = 10
# N_end = 100
# dt = 0.05

# tt = 0:dt:N_end
# fig = Figure()
# fig2 = Figure()
# fig3 = Figure()
# ax = Axis( fig[1, 1] )
# ax2 = Axis( fig2[1, 1] )
# ax3 = Axis( fig3[1, 1] )
# for _ in 1:20
#     local xx = OU_solve( theta, x0, N_end; dt = dt )
#     lines!( ax, tt, xx[1, :] )
#     diffs = calculate_diffs( xx, [1, 2], dt )
#     lines!( ax2, tt[1:end-2], diffs[1][1, :] )
#     lines!( ax3, tt[1:end-4], diffs[2][1, :] )
# end

# display(fig)
# display(fig2)
# display(fig3)

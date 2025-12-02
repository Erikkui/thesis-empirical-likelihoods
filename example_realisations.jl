include( "OU/OU_solve.jl" )
include( "L3/lorenz3.jl" )
include( "blowfly/blowfly_solve.jl" )

using CairoMakie
using DifferentialEquations
using Distributions
colors = Makie.to_colormap(:seaborn_dark)
set_theme!(palette = (; color = colors))

## OU
theta = [3, 2]
x0 = 4
Ndata = 200
dt = 0.01

tt = ( 0:Ndata ) .* dt
fig = Figure( size = (600, 225), padding = 1, fontsize = 13 )
ax = Axis( fig[1, 1], height = 200 )
for _ in 1:3
    local xx = OU_solve( theta, x0, Ndata; dt = dt )
    lines!( ax, tt, xx[1, :], linewidth=2.5 )
end
# ax.xticklabelsize = 20
# ax.yticklabelsize = 20

save("OU_simulation.pdf", fig)
display(fig)

# L3
# theta = [10.0, 28.0, 8/3]
# init = [12.577, 19.471, 23.073]
# dt = 0.01;
# Ndata = 5000 + 10;

# R = lorenz_solve( init, theta, Ndata; dt = dt )
# t = ( 0:Ndata ) .* dt
# t = t[11:end]

# fig2 = Figure( size = (700.33, 300), fontsize = 16, figure_padding = (2, 21, 20, 0) )
# # fig3 = Figure( size = (266.67, 266.67), fontsize = 16, figure_padding = 1 )
# ax_phase = Axis( fig2[1, 1:6], height = 220 )
# ax_3d = Axis3( fig2[1, 7:10], azimuth = -8*pi/19, xlabeloffset = 25, zlabeloffset = 35, height = 220)
# lines!( ax_3d, R[1, :], R[2, :], R[3, :], linewidth=1 )
# tind = 1:1500
# labels = ["x", "y", "z"]
# for ii in 1:3
#     lines!( ax_phase, t[tind], R[ii, tind], linewidth=2, label = labels[ii] )
# end
# # axislegend( ax_phase; position = :cb, bgcolor = (:white, 0.8), framevisible = false )
# axislegend( ax_phase; position = (0.5, 0.025), framevisible = false )

# display(fig2)
# # display(fig3)
# save("L3_phase_simulation.pdf", fig2)
# save("L3_3d_simulation.pdf", fig3)

# Blowfly model
# t_max = 400
# tt = range( 0, t_max, t_max )  # Time vector from 0 to N with N+1 points

# delta = 0.16
# P = 6.5
# N_0 = 400
# sigma_p = 0.1
# tau = 14
# sigma_d = 0.1
# theta = [ delta, P, N_0, sigma_p, tau, sigma_d ]

# init = 180
# burn_in = 0


# fig = Figure( size = (600, 250), padding = 1, fontsize = 13 )

# ax = Axis( fig[1, 1],
#     # xgridvisible = true,
#     # ygridvisible = true,
#     xticks = 0:50:t_max,
#     yticks = 0:2000:6000
# )

# for nn = 1:3
#     local N_sim, _ = blowfly_solve( theta, t_max, N_init = init, burn_in = burn_in )
#     local tt = collect( 0:size( N_sim, 2) - 1 )
#     lines!( ax, tt[:], N_sim[:], linewidth = 2.5)
# end

# #

# save( "blowfly_simulation.pdf", fig )
# display( fig )

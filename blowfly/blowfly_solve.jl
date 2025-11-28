
function blowfly_solve( theta, t; N_init = 180, burn_in = 20, mu = 1.0 )
    # Unpack parameters
    delta, P, N_0, sigma2_p, tau, sigma2_d = theta

    # Initializing needed variables
    lag = Int( round(tau) ) > 0 ? Int( round(tau) ) : 1     # Lag time in days
    total_time = t + lag + burn_in
    N = Matrix{Float64}(undef, 1, total_time )  # Population size vector
    N[1, 1:lag] .= N_init  # Initializing first lag days with N_init

    gamma_p = Gamma( mu^2/sigma2_p, sigma2_p/mu )   # Gamma distribution for reproduction rate
    gamma_d = Gamma( mu^2/sigma2_d, sigma2_d/mu )   # Gamma distribution for death rate
    ee = rand( gamma_p, total_time )                # Reproduction rate noise
    epsilon = rand( gamma_d, total_time )           # Death rate noise

    # Iterate over time steps; first burn in period, then actual simulation
    # Note: We start from lag+1 because we need to access N[ii - lag]
    for ii in lag+1:total_time
        Nlag = @view N[ 1, ii - lag ]       # Population size at lag time
        Nprev = @view N[ 1, ii - 1 ]        # Population size at previous time step
        ee_t = @view ee[ ii - 1 ]           # Reproduction rate noise at lag time
        epsilon_t = @view epsilon[ ii - 1 ] # Death rate noise at lag time

        # Expected values: R ~ Poisson, S ~ binom
        R_t = P*Nlag*exp.( -Nlag/N_0 ).*ee_t      # E(X) = lambda
        S_t = Nprev*exp.( -delta*epsilon_t )      # E(X) = n*p
        N[ii] = R_t + S_t[]

    end

    # Remove burn-in period and lag padding
    N_burned = N[ 1:1, 1+lag+burn_in:end ]

    return N, N_burned
end


## For testing the blowfly model
# using Distributions
# using CairoMakie
# N = 500
# # tt = range( 0, N, N )  # Time vector from 0 to N with N+1 points

# delta = 0.16
# P = 6.5
# N_0 = 400
# sigma_p = 0.1
# tau = 14
# sigma_d = 0.1
# theta = [ delta, P, N_0, sigma_p, tau, sigma_d ]

# init = 180
# burn_in = 50

# N_sim, N_burned = blowfly_solve( theta, N, N_init = init, burn_in = burn_in )

# println( typeof(N_sim), " ", size(N_sim), " ", typeof(N_burned), " ", size(N_burned) )

# fig = Figure( size = (1200, 600) )

# ax = Axis( fig[1, 1],
#     title = "Blowfly Population Dynamics",
#     xlabel = "Time (days)",
#     ylabel = "Population Size",
#     xgridvisible = true,
#     ygridvisible = true,
#     xticks = 0:50:N,
# )

# for nn = 1:3
#     N_sim, _ = blowfly_solve( theta, N, N_init = init, burn_in = burn_in )
#     tt = collect(0:size( N_sim, 2) - 1)
#     lines!( ax, tt[1:tau+burn_in], N_sim[1, 1:tau+burn_in], linewidth = 1.5, linestyle=:dash )
#     lines!( ax, tt[1+tau+burn_in:end], N_sim[1, 1+tau+burn_in:end], linewidth = 1.5 )

# end

# # 554
# theta2 = [0.1876    8.0774  400    0.2900   19.7181    0.2943]

# N_sim, _ = blowfly_solve( theta2, N, N_init = init, burn_in = burn_in )
# N_len = size( N_sim, 2 )

# lines!(ax, 0:N_len-1, N_sim[1, :], linewidth = 3, color=:red )


# display( fig )

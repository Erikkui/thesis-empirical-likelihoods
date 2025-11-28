function create_synth_data( data::Dict{Symbol, Any} )
    # Create synthetic data for the blowfly model
    use_diff = data[:use_diff]
    use_2D = data[:use_2D]
    diff_order = data[:diff_order]

    params = copy(data[:params])
    init = data[:synth_init]
    dt = data[:synth_dt]
    N_end = data[:N_end]
    nepo = data[:nepo]

    R0 = Vector{ Matrix{Float64} }(undef, nepo)

    if use_diff == 1 && ~use_2D
        R0_diff = Vector{ Vector{Matrix{Float64}} }(undef, length(diff_order) )
        for ii in eachindex(diff_order)
            R0_diff[ii] = Vector{ Matrix{ Float64 } }(undef, 0)
        end
    end

    for ii in 1:nepo
        N = lorenz_solve( init, params, N_end; dt = dt )
        R0[ii] = N
        if use_diff == 1
            R0_diff_ii = calculate_diffs( N, diff_order, dt, use_2D = use_2D )
            if use_2D
                R0[ii] = vcat( N, R0_diff_ii... )  # Store the population size with differences
            else
                push!.( R0_diff, R0_diff_ii )
            end
        end
    end

    data[:nobs] = size( R0[1], 2 )
    data[:R0] = R0

    if use_diff == 1 && use_2D != 1
        data[:R0_diff] = R0_diff
    end

    return data
end






# Version when using 2D arrays for R0 and R0_diff
# function create_synth_data_2D( data::Dict{Symbol, Any} )
#     # Create synthetic data for the blowfly model
#     use_diff = data[:use_diff]
#     params = copy(data[:theta])
#     init = data[:synth_init]
#     burn_in = data[:synth_burn_in]
#     t = data[:synth_t]
#     nepo = data[:nepo]

#     R0 = Vector{ Matrix{Float64} }(undef, nepo)
#     R0_full = Vector{ Matrix{Float64} }(undef, nepo)

#     if use_diff == 1
#         # If using differences, we need to adjust the time vector
#         # t = range( t[1], t[end]+1, length(t)+1 )
#         data[:t_diff] = t
#         R0_diff = Vector{ Matrix{Float64} }(undef, nepo)
#     end

#     for ii in 1:nepo
#         N, N_burned = blowfly_solve( params, t, N_init = init, burn_in = burn_in )
#         R0_full[ii] = N
#         if use_diff == 1
#             N_diff_prev = N_burned[1:1, 1:end-2]  # Calculating central difference
#             N_diff_next = N_burned[1:1, 3:end]
#             N_diff = (N_diff_next - N_diff_prev) ./ (2*t.step)  # Central difference
#             N_burned = N_burned[1:1, 2:end-1]
#             # R0_diff[ii] = N_diff
#         end
#         R0[ii] = [ N_burned; N_diff ]  # Store the population size without the last time point
#     end

#     data[:nobs] = size( R0[1], 2 )
#     data[:R0] = R0
#     data[:R0_full] = R0_full

#     if use_diff == 1
#         data[:R0_diff] = R0_diff
#     end

#     return data
# end

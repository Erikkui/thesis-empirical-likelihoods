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
        N = OU_solve( params, init, N_end; dt = dt )
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

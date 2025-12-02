function create_synth_data( data::Dict{Symbol, Any} )
    # Create synthetic data for the blowfly model
    use_diff = data[:use_diff]
    diff_order = data[:diff_order]

    params = copy(data[:params])
    init = data[:synth_init]
    Ndata = data[:Ndata]
    dt = data[:synth_dt]
    nepo = data[:nepo]

    R0 = Vector{ Matrix{Float64} }(undef, nepo)

    if use_diff == 1
        R0_diff = Vector{ Vector{Matrix{Float64}} }(undef, length(diff_order) )
        for ii in eachindex(diff_order)
            R0_diff[ii] = Vector{ Matrix{ Float64 } }(undef, 0)
        end
    end

    for ii in 1:nepo
        N = OU_solve( params, init, Ndata; dt = dt )
        R0[ii] = N
        if use_diff == 1
            R0_diff_ii = calculate_diffs( N, diff_order, dt)
            push!.( R0_diff, R0_diff_ii )
        end
    end

    data[:nobs] = size( R0[1], 2 )
    data[:R0] = R0

    if use_diff == 1
        data[:R0_diff] = R0_diff
    end

    return data
end

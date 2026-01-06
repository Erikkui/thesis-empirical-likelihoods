# For cases where difference and no 2D data is used
function Wrapper(theta::Vector{Float64}, data::Dict{Symbol, Any})
    dt = data[:synth_dt]
    t = data[:synth_N]
    burn_in = data[:synth_burn_in]
    init = data[:synth_init]
    R0_all = data[:R0_all]
    case_nsim, _ = data[:case_dim]
    chamfer = data[:chamfer]
    chamfer_k = data[:chamfer_k]
    eCDF = data[:eCDF]
    bins = data[:bins]
    nL = data[:nL]
    diff_order = data[:diff_order]
    use_diff = data[:use_diff]

    # early return if delay parameter is definitely too  large (blowfly lifetime max one month)
    if theta[5] > 60
        return data, nothing
    end

    ######## Initialize data structures; as undef first for conactenation to work
    chamfer_dists = Array{Float64}(undef, case_nsim, 0)
    cdfs = Array{Float64}(undef, case_nsim, 0)
    if chamfer == 1
        chamfer_dists = zeros( case_nsim, length(chamfer_k)*( use_diff*length(diff_order) + 1 ) )
    end
    if eCDF == 1
        cdfs = zeros( case_nsim, ( 1 + use_diff*length(diff_order) )*data[:nbin]*nL )
    end
# println( chamfer_dists)
    ######### Simulations for proposal parameter
    R_sim_all = Vector{ Matrix{Float64} }(undef, 0 )
    if use_diff == 1
        R_sim_diff_all = Vector{ Vector{ Matrix{Float64} } }(undef, length( diff_order ) )
        for ii in eachindex( diff_order )
            R_sim_diff_all[ii] = Vector{ Matrix{Float64} }(undef, 0 )
        end
    end

    # Generate simulated data
    for _ in 1:case_nsim
        N, _ = blowfly_solve( theta, t, N_init = init, burn_in = burn_in )
        push!( R_sim_all, N )
        if use_diff == 1
            N_diffs = calculate_diffs( N, diff_order, dt )
            for ii in eachindex( diff_order )
                push!( R_sim_diff_all[ii], N_diffs[ii] )
            end
        end
    end

    if use_diff == 1
        R_sim_all = [ R_sim_all, R_sim_diff_all... ]  # Concatenate all simulated data into vector of matrices
    else
        R_sim_all = [ R_sim_all ]
    end
# println( chamfer_dists)
    # Create eCDFs and chamfer distances for each simulations
    for dd in eachindex(R0_all) # Loop over data types: original and differences
        x = R0_all[dd]
        R_sim_dd = R_sim_all[dd]
        data[:bins] = bins[dd]

        for ii in eachindex( R_sim_dd )
            y = copy( R_sim_dd[ii] )
            data, cdfs_ii, chamfer_ii = resample_data( x, y, data )

            if size( cdfs_ii, 1 ) > 1
                cdfs_ii = mean( cdfs_ii, dims=1 )
                chamfer_ii = mean( chamfer_ii, dims=1 )
            end

            if eCDF == 1
                cdfs_inds = (dd-1)*data[:nbin]*nL+1:dd*data[:nbin]*nL
                cdfs[ ii, cdfs_inds ] = cdfs_ii
            end
            if chamfer == 1
                chamf_inds = (dd-1)*length( chamfer_k )+1:dd*length( chamfer_k )
                chamfer_dists[ ii, chamf_inds ] = chamfer_ii
            end
        end
    end
    data[:bins] = bins  # Reset bins in data dict

    # Normalize chamfer
    means = mean(chamfer_dists, dims = 1)
    stds = std( chamfer_dists, dims = 1)
    chamfer_dists = (chamfer_dists .- means) ./ stds

    # println( chamfer_dists)
    summary_stats = [ cdfs chamfer_dists ]

    if data[:log] == "log"
        summary_stats = log.(summary_stats)
    end

    muu = mean( summary_stats, dims=1 )

    # Compute muu and C based on simulated data
    if data[:C_how] == "don"
        C = donsker( muu, data[:nobs] )
    elseif data[:C_how] == "cov"
        C = cov( summary_stats )
    end

    data[:C] = C

    return data, muu
end

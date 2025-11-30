# For cases where difference and no 2D data is used
function Wrapper(theta::Vector{Float64}, data::Dict{Symbol, Any})
    dt = data[:synth_dt]
    t = data[:synth_N]
    burn_in = data[:synth_burn_in]
    init = data[:synth_init]
    R0_all = data[:R0_all]
    nepo = data[:nepo]
    nobs = data[:nobs]
    case_nsim, case_nrep = data[:case_dim]
    chamfer = data[:chamfer]
    chamfer_k = data[:chamfer_k]
    eCDF = data[:eCDF]
    LL = data[:LL]
    nL = data[:nL]
    create_kdtree = data[:create_kdtree]
    kn = data[:kn]
    case = data[:case]
    data_dim = data[:data_dim]
    diff_order = data[:diff_order]
    use_diff = data[:use_diff]



    ######## Initialize data structures; as undef first for conactenation to work
    chamfer_dists = Array{Float64}(undef, case_nsim, 0)
    cdfs = Array{Float64}(undef, case_nsim, 0)
    if chamfer == 1
        chamfer_dists = zeros( case_nsim, length(chamfer_k)*( use_diff*length(diff_order) + 1 ) )
    end
    if eCDF == 1
        cdfs = zeros( case_nsim, ( 1 + use_diff*length(diff_order) )*data[:nbin]*nL )
    end


    ######### Simulations for proposal parameter
    R_sim_all = Matrix{Float64}(undef, case_nsim, nobs)
    if use_diff == 1
        R_sim_diff_all = Vector{ Matrix{Float64} }(undef, length( diff_order ) )
        for ii in eachindex( diff_order )
            R_sim_diff_all[ii] = Matrix{ Float64 }(undef, case_nsim, nobs - 2*diff_order[ii] )
        end
    end

    for sim in 1:case_nsim
        N, _ = blowfly_solve( theta, t, N_init = init, burn_in = burn_in )
        R_sim_all[ sim, : ] = N
        if use_diff == 1
            N_diffs = calculate_diffs( N, diff_order, dt )
            for ii in eachindex( diff_order )
                R_sim_diff_all[ii][ sim, : ] = N_diffs[ii]
            end
        end
    end

    if use_diff == 1
        R_sim_all = [ R_sim_all, R_sim_diff_all... ]  # Concatenate all simulated data into vector of matrices
    else
        R_sim_all = [ R_sim_all ]
    end

    # Create eCDFs and chamfer distances for each simulations
    for dd in eachindex(R0_all)
        x = R0_all[dd]
        R_sim_dd = R_sim_all[dd]

        for (ii, y) in enumerate( eachrow( R_sim_dd ) )
            if -1 in data[:LL]          # Signal feature
                cdfs_ii, chamfer_ii = create_summaries( y, data )
                LL = LL[ 2:end ]  # Remove -1 for next features
                cdfs_ii_2, chamfer_ii_2 = create_summaries( x, y, data, create_kdtree, kn )
                LL = data[:LL]  # Reset LL

                cdfs_ii = hcat( cdfs_ii, cdfs_ii_2 )
                chamfer_ii = hcat( chamfer_ii, chamfer_ii_2 )


            else
                cdfs_ii, chamfer_ii = create_summaries( x, y, data, create_kdtree, kn )
            end

            if eCDF == 1
                cdfs[ ii, : ] = cdfs_ii
            end
            if chamfer == 1
                chamfer_dists[ ii, : ] = chamfer_ii
            end

        end
    end



    summary_stats = [ cdfs chamfer_dists ]

    if data[:log] == "log"
        summary_stats = log.(summary_stats)
    end

    # Compute muu and C based on simulated data
    if data[:case] == "bsl"

        if data[:C_how] == "don"
            muu = mean( summary_stats, dims=1 )
            C = donsker( muu, data[:nobs] ) #TODO check nobs
        elseif data[:C_how] == "cov"
            C = cov( summary_stats )
        end

        data[:C] = C
    end

    if size( summary_stats, 1 ) > 1
        summary_stats = mean( summary_stats, dims=1 )
    end

    return data, summary_stats
end

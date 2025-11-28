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
    nL = data[:nL]
    create_kdtree = data[:create_kdtree]
    kn = data[:kn]
    case = data[:case]
    data_dim = data[:data_dim]
    diff_order = data[:diff_order]
    use_diff = data[:use_diff]
    signal_only = data[:signal_only]

    # If using only signal and no dists: make ecdf from each simulated signal
    if signal_only
        case_nrep = case_nsim   # Update nrep for array initialization and indices calc
    end

    ######## Initialize data structures; as undef first for conactenation to work
    chamfer_dists = Array{Float64}(undef, case_nrep, 0)
    cdfs = Array{Float64}(undef, case_nrep, 0)
    if chamfer == 1
        chamfer_dists = zeros( case_nrep, length(chamfer_k)*( use_diff*length(diff_order) + 1 ) )
    end
    if eCDF == 1
        if signal_only && case == "bsl"
            cdfs = zeros( nepo, ( 1 + use_diff*length(diff_order) )*data[:nbin]*nL )
            chamfer_dists = Array{Float64}(undef, nepo, 0)    # For concatenation to work correctly
        else
            cdfs = zeros( case_nrep, ( 1 + use_diff*length(diff_order) )*data[:nbin]*nL )
        end
    end

    # For simulation based muu and C
    chamfer_dists_hat = Array{Float64}(undef, case_nrep, 0)
    cdfs_hat = Array{Float64}(undef, case_nrep, 0)
    if case == "bsl"
        if chamfer == 1
            chamfer_dists_hat = zeros( case_nrep, length(chamfer_k)*( use_diff*length(diff_order) + 1 ) )
        end
        if eCDF == 1
            cdfs_hat = zeros( case_nrep, ( 1 + use_diff*length(diff_order) )*data[:nbin]*nL )
        end
    else    # Save bins for GSL as they do not change.
        binss_R = copy( data[:binss] )
    end


    ######### Simulations for proposal parameter
    R_sim_all = Matrix{Float64}(undef, data_dim, nobs*case_nsim)
    R_sim_diff_all = Matrix{Float64}(undef, data_dim, 0)
    if use_diff == 1
        R_sim_diff_all = [ Matrix{Float64}(undef, data_dim, (nobs - 2*ii)*case_nsim) for ii in 1:length(diff_order) ]
    end
    for sim in 1:case_nsim
        ind_start = (sim - 1)*nobs + 1
        ind_end = sim*nobs
        _, N_burned = blowfly_solve( theta, t, N_init = init, burn_in = burn_in )
        R_sim_all[ :, ind_start:ind_end ] = @view N_burned[ 1:data_dim, : ]

        if use_diff == 1
            N_diffs = calculate_diffs( N_burned, diff_order, dt )
            for ii in eachindex( diff_order )
                order = diff_order[ii]
                ind_start = (sim - 1)*(nobs - 2*order) + 1
                ind_end = sim*(nobs - 2*order)
                R_sim_diff_all[ii][ 1, ind_start:ind_end ] = @view N_diffs[ii][ 1, : ]
            end
        end
    end

    if use_diff == 1
        R_sim_all = [ R_sim_all, R_sim_diff_all... ]  # Concatenate all simulated data into vector of matrices
    else
        R_sim_all = [ R_sim_all ]
    end

    # Signal only as a special case
    if signal_only
        cdf_ii = Array{Float64}(undef, 1, 0)
        for dd in eachindex( R0_all )  # Loop through each data set: signal and its difference(s)
            Rsim = R_sim_all[dd]
            R = R0_all[dd]

            ndata_Rsim = Int( size( Rsim, 2 ) / case_nsim )
            ndata_R = Int( size( R, 2 ) / nepo )
            cdfs_inds = (dd-1)*data[:nbin]*nL+1:dd*data[:nbin]*nL

            if case == "gsl"
                # Pick bins for data/diff data
                binss_inds = (dd-1)*nL+1:dd*nL
                data[:binss] = binss_R[ binss_inds ]

                for rep in 1:case_nrep
                    ind_start = (rep-1)*ndata_Rsim+1
                    ind_end = rep*ndata_Rsim
                    Y = @view Rsim[ :, ind_start:ind_end ]

                    cdf_ii, _ = create_summaries( Y, data )
                    cdfs[ rep, cdfs_inds ] = cdf_ii
                end
            else

                # Create bins based on simulated data
                dataL = copy( data )
                r = Array{Float64}(undef, size( Rsim, 1 ), nepo*ndata_R )
                for kk=1:nepo
                    r[ :, (kk-1)*ndata_Rsim+1:kk*ndata_Rsim ] = @view Rsim[ :, (kk-1)*ndata_Rsim+1:kk*ndata_Rsim ]
                end
                dataL[:r] = r
                m = maximum( minimum( Rsim, dims=2 ) )
                M = minimum( maximum( Rsim, dims=2 ) )
                dataL[:minmax] = [ m, M ]
                data[:binss] = [ bin_select( dataL ) ]

                for ii in 1:nepo
                    Rii = @view R[ :, (ii-1)*ndata_R+1:ii*ndata_R ]
                    cdf_ii, _ = create_summaries( Rii, data )  # Summaries from data
                    cdfs[ ii, cdfs_inds ] = cdf_ii
                end

                for rep in 1:case_nrep
                    ind_start = (rep-1)*ndata_Rsim+1
                    ind_end = rep*ndata_Rsim
                    Ysim = @view Rsim[ :, ind_start:ind_end ]

                    cdf_ii_hat, _ = create_summaries( Ysim, data ) # Summaries from sim data; used to create mu and C
                    cdfs_hat[ rep, cdfs_inds ] = cdf_ii_hat
                end

            end
        end

    else

        # For each data set: signal and its difference(s)
        for dd in eachindex(R0_all)
            R0 = R0_all[dd]
            Rsim = R_sim_all[dd]

            ndata_R0 = size( R0, 2 )
            ndata_Rsim = size( Rsim, 2 )

            if case == "bsl"
                if eCDF == 1
                    data[:binss] = Vector{Vector{Float64}}(undef, 1)
                    data, _ = resample_data( Rsim, data )  # Update bins according to simulated data
                end
            end

            for rep in 1:case_nrep

                # Randomly sampling from the trajectories
                if case == "gsl"
                    x_ind = rand( 1:ndata_R0, ndata_R0 )
                    unique!( x_ind )
                    x = @view R0[ :, x_ind ]
                    y = @view Rsim[:, :]

                    # Pick bins for data/diff data
                    binss_inds = (dd-1)*nL+1:dd*nL
                    data[:binss] = @view binss_R[ binss_inds ]

                else

                    x_ind = rand( 1:ndata_Rsim, ndata_Rsim )
                    unique!( x_ind )
                    x = @view Rsim[ :, x_ind ]
                    y = @view R0[:, :]

                    # For muu based on data
                    yhat_inds = setdiff( 1:ndata_Rsim, x_ind )
                    yhat = @view Rsim[ :, yhat_inds ]

                end

                cdf_ii, c_ii = create_summaries( x, y, data, create_kdtree, kn )
                if chamfer == 1
                    chamf_inds = (dd-1)*length( chamfer_k )+1:dd*length( chamfer_k )
                    chamfer_dists[ rep, chamf_inds ] = c_ii
                end
                if eCDF == 1
                    cdfs_inds = (dd-1)*data[:nbin]*nL+1:dd*data[:nbin]*nL
                    cdfs[ rep, cdfs_inds ] = cdf_ii
                end

                # Calculate summaries for muu and C based on simulated data
                if case == "bsl"
                    cdf_ii, c_ii = create_summaries( x, yhat, data, create_kdtree, kn )
                    if chamfer == 1
                        chamf_inds = (dd-1)*length( chamfer_k )+1:dd*length( chamfer_k )
                        chamfer_dists_hat[ rep, chamf_inds ] = c_ii
                    end
                    if eCDF == 1
                        cdfs_inds = (dd-1)*data[:nbin]*nL+1:dd*data[:nbin]*nL
                        cdfs_hat[ rep, cdfs_inds ] = cdf_ii
                    end

                    data[:binss] = Vector{Vector{Float64}}(undef, nL)
                end

            end
        end
    end

    summary_stats = [ cdfs chamfer_dists ]

    if data[:log] == "log"
        summary_stats = log.(summary_stats)
    end

    if size( summary_stats, 1 ) > 1
        summary_stats = mean( summary_stats, dims=1 )
    end

    # Compute new muu and C based on simulated data
    if data[:case] == "bsl"
        summary_stats_hat = [ cdfs_hat chamfer_dists_hat ]

        if data[:log] == "log"
            summary_stats_hat = log.(summary_stats_hat)
        end

        muu = mean( summary_stats_hat, dims=1 )

        if data[:C_how] == "don"
            C = donsker( muu, data[:N] )
        elseif data[:C_how] == "cov"
            C = cov( summary_stats_hat )
        end

        data[:muu] = muu
        data[:C] = C

        # println("Summary Stats Hat: ", summary_stats_hat)
        # println("C: ", C, "\n")

    else
        data[:binss] = binss_R  # Put original bins back for next round
    end

    return data, summary_stats
end

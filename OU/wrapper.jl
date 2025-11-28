# For cases where difference and no 2D data is used
function Wrapper(theta::Vector{Float64}, data::Dict{Symbol, Any})
    dt = data[:synth_dt]
    N_end = data[:N_end]
    init = data[:synth_init]
    R0_all = data[:R0_all]
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
        cdfs = zeros( case_nrep, ( 1 + use_diff*length(diff_order) )*data[:nbin]*nL )
    end

    # For simulation based muu and C
    chamfer_dists_hat = Array{Float64}(undef, case_nrep, 0)
    cdfs_hat = Array{Float64}(undef, case_nrep, 0)
    if case == "bsl"
        if chamfer == 1
            chamfer_dists_hat = similar( chamfer_dists )
        end
        if eCDF == 1
            cdfs_hat = similar( cdfs )
        end
    else    # Save bins for GSL as they do not change.
        binss_R = data[:binss]
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

        N = OU_solve( theta, init, N_end; dt = dt )
        R_sim_all[ :, ind_start:ind_end ] = N

        if use_diff == 1
            N_diffs = calculate_diffs( N, diff_order, dt )
            for ii in eachindex( diff_order )
                order = diff_order[ii]
                ind_start = (sim - 1)*(nobs - 2*order) + 1
                ind_end = sim*(nobs - 2*order)
                R_sim_diff_all[ii][ 1, ind_start:ind_end ] = @view N_diffs[ii][ 1, : ]
            end
        end
    end


    if use_diff == 1
        R_sim_all = [ R_sim_all, R_sim_diff_all... ]  # Concatenate all simulated data
    else
        R_sim_all = [ R_sim_all ]
    end


    # Signal only as a special case
    if signal_only
        if case == "gsl"
            for dd in eachindex(R_sim_all)
                Rsim = R_sim_all[dd]
                ndata_Rsim = Int( size( Rsim, 2 ) / case_nsim )

                # Pick bins for data/diff data
                binss_inds = (dd-1)*nL+1:dd*nL
                data[:binss] = binss_R[ binss_inds ]

                # Both x and y needed for functions to work

                for rep in 1:case_nrep
                    ind_start = (rep-1)*ndata_Rsim+1
                    ind_end = rep*ndata_Rsim
                    Y = @view Rsim[ :, ind_start:ind_end-2 ]

                    cdf_ii, _ = create_summaries( Y, data )

                    cdfs_inds = (dd-1)*data[:nbin]*nL+1:dd*data[:nbin]*nL
                    cdfs[ rep, cdfs_inds ] = cdf_ii
                end
            end
        else
            # TODO BSL case when using only signal; ask Heikki for details
        end

    else

        # For each data set: signal and its difference(s)
        for dd in eachindex(R0_all)
            R0 = R0_all[dd]
            Rsim = R_sim_all[dd]

            ndata_R0 = size( R0, 2 )
            ndata_Rsim = size( Rsim, 2 )

            for rep in 1:case_nrep

                # Randomly sampling from the trajectories
                if case == "gsl"
                    x_ind = rand( 1:ndata_R0, ndata_R0 )
                    unique!( x_ind )
                    x = @view R0[ :, x_ind ]
                    y = @view Rsim[:, :]

                    # Pick bins for data/diff data
                    binss_inds = (dd-1)*nL+1:dd*nL
                    data[:binss] = binss_R[ binss_inds ]
                else
                    if eCDF == 1
                        data[:binss] = Vector{Vector{Float64}}(undef, 1)
                        data, _ = resample_data( Rsim, data )  # Update bins for each rep
                    end
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

    else
        data[:binss] = binss_R  # Put original bins back for next round
    end

    return data, summary_stats
end


# For cases where no difference is used
# function Wrapper(theta, data )
#     t = data[:synth_t]
#     init = data[:synth_init]
#     burn_in = data[:synth_burn_in]
#     R0_all = data[:R0_all]
#     nobs = data[:nobs]
#     nepo = data[:nepo]
#     case_nsim, case_nrep = data[:case_dim]
#     chamfer = data[:chamfer]
#     eCDF = data[:eCDF]
#     LL = data[:LL]
#     nL = data[:nL]
#     create_kdtree = data[:create_kdtree]
#     kn = data[:kn]


#     ######## Initialize data structures; as undef first for conactenation to work
#     chamfer_dists = Array{Float64}(undef, case_nrep, 0)
#     cdfs = Array{Float64}(undef, case_nrep, 0)
#     if chamfer == 1
#         chamfer_k = data[:chamfer_k]
#         chamfer_dists = zeros( case_nrep, length( chamfer_k ) )
#     end
#     if eCDF == 1
#         cdfs = zeros( case_nrep, data[:nbin]*nL )
#     end


#     # For simulation based muu and C
#     chamfer_dists_hat = Array{Float64}(undef, case_nrep, 0)
#     cdfs_hat = Array{Float64}(undef, case_nrep, 0)
#     if data[:case] == "bsl"
#         if chamfer == 1
#             chamfer_dists_hat = similar( chamfer_dists )
#         end
#         if eCDF == 1
#             cdfs_hat = similar( cdfs )
#         end
#     end


#     ######### Simulations for proposal parameter
#     R_sim_all = Matrix{Float64}(undef, length(init), nobs*case_nsim)
#     for sim in 1:case_nsim
#         ind_start = (sim - 1)*nobs + 1
#         ind_end = sim*nobs
#         _, N_burned = blowfly_solve( theta, t, N_init = init, burn_in = burn_in )
#         # println( size( N_burned ) )

#         R_sim_all[ :, ind_start:ind_end ] = N_burned
#     end


#     if data[:resample] == 1

#         for rep in 1:case_nrep

#             # Randomly sampling from the trajectories
#             if data[:case] == "gsl"
#                 x_ind = rand( 1:nepo*nobs, nepo*nobs)
#                 unique!( x_ind )
#                 x = @view R0_all[ :, x_ind ]
#                 y = R_sim_all
#             else
#                 data, _ = resample_data( R_sim_all, data )  # Update bins for each rep
#                 x_ind = rand( 1:case_nsim*nobs, case_nsim*nobs )
#                 unique!( x_ind )
#                 x = @view R_sim_all[ :, x_ind ]
#                 y = R0_all

#                 # For muu based on data
#                 yhat_inds = setdiff( 1:case_nsim*nobs, x_ind )
#                 yhat = @view R_sim_all[ :, yhat_inds ]
#             end

#             cdf_ii, c_ii = create_summaries( x, y, data, create_kdtree, kn )
#             cdfs[ rep, : ] = cdf_ii
#             chamfer_dists[ rep, : ] = c_ii

#             # Calculate summaries for muu and C based on simulated data
#             if data[:case] == "bsl"
#                 cdf_ii, c_ii = create_summaries( x, yhat, data, create_kdtree, kn )
#                 cdfs_hat[ rep, : ] = cdf_ii
#                 chamfer_dists_hat[ rep, : ] = c_ii
#                 data[:binss] = Vector{Vector{Float64}}(undef, nL)
#             end

#         end

#     else
#         # Resample by epochs
#         # R = Vector{Any}(undef, 2)
#         # R[1] = y

#         # for rep in 1:case_nrep
#         #     row_ind = (samp-1)*case_nrep + rep
#         #     idata = rand(1:nepo)
#         #     R[2] = data[:R0][idata]
#         #     data, cdf_row = pairwise_distL(R, data)
#         #     cdfs[row_ind, :] = cdf_row
#         # end

#     end


#     summary_stats = [ cdfs chamfer_dists ]

#     if data[:log] == "log"
#         summary_stats = log.(summary_stats)
#     end

#     if size( summary_stats, 1 ) > 1
#         summary_stats = mean( summary_stats, dims=1 )
#     end

#     # Compute new muu and C based on simulated data
#     if data[:case] == "bsl"
#         summary_stats_hat = [ cdfs_hat chamfer_dists_hat ]

#         if data[:log] == "log"
#             summary_stats_hat = log.(summary_stats_hat)
#         end

#         muu = mean( summary_stats_hat, dims=1 )

#         if data[:C_how] == "don"
#             C = donsker( muu, data[:N] )
#         elseif data[:C_how] == "cov"
#             C = cov( summary_stats_hat )
#         end

#         data[:muu] = muu
#         data[:C] = C
#     end

#     return data, summary_stats
# end


# Wrapper fot testing 2D arrays
function Wrapper2D(theta::Vector{Float64}, data::Dict{Symbol, Any}, use_diff::Int64)
    t = data[:synth_t]
    burn_in = data[:synth_burn_in]
    init = data[:synth_init]
    R0_all = data[:R0_all]
    nobs = data[:nobs]
    nepo = data[:nepo]
    case_nsim, case_nrep = data[:case_dim]
    chamfer = data[:chamfer]
    chamfer_k = data[:chamfer_k]
    eCDF = data[:eCDF]
    LL = data[:LL]
    nL = data[:nL]
    create_kdtree = data[:create_kdtree]
    kn = data[:kn]
    case = data[:case]

    R0_dim, R0_ndata = size( R0_all )

    ######## Initialize data structures; as undef first for conactenation to work
    chamfer_dists = Array{Float64}(undef, case_nrep, 0)
    cdfs = Array{Float64}(undef, case_nrep, 0)
    if chamfer == 1
        chamfer_dists = zeros( case_nrep, length( chamfer_k ) )
    end
    if eCDF == 1
        cdfs = zeros( case_nrep, data[:nbin]*nL )
    end


    # For simulation based muu and C
    chamfer_dists_hat = Array{Float64}(undef, case_nrep, 0)
    cdfs_hat = Array{Float64}(undef, case_nrep, 0)
    if case == "bsl"
        if chamfer == 1
            chamfer_dists_hat = similar( chamfer_dists )
        end
        if eCDF == 1
            cdfs_hat = similar( cdfs )
        end
    end


    ######### Simulations for proposal parameter
    R_sim_all = Matrix{Float64}(undef, R0_dim, nobs*case_nsim)
    for sim in 1:case_nsim
        ind_start = (sim - 1)*nobs + 1
        ind_end = sim*nobs

        _, N_burned = blowfly_solve( theta, t, N_init = init, burn_in = burn_in )
        N_diff = calculate_diffs( N_burned, data[:diff_order], data[:use_2D] )

        R_sim_all[ :, ind_start:ind_end ] = vcat( N_burned, N_diff... )
    end

    ndata_Rsim = size( R_sim_all, 2 )

    for rep in 1:case_nrep

        # Randomly sampling from the trajectories
        if case == "gsl"
            x_ind = rand( 1:R0_ndata, R0_ndata )
            unique!( x_ind )
            x = @view R0_all[ :, x_ind ]
            y = R_sim_all
        else
            data, _ = resample_data( R_sim_all, data )  # Update bins for each rep
            x_ind = rand( 1:ndata_Rsim, ndata_Rsim )
            unique!( x_ind )
            x = @view R_sim_all[ :, x_ind ]
            y =  R0_all

            # For muu based on data
            yhat_inds = setdiff( 1:ndata_Rsim, x_ind )
            yhat = @view R_sim_all[ :, yhat_inds ]
        end

        cdf_ii, c_ii = create_summaries( x, y, data, create_kdtree, kn )
        if chamfer == 1
            # chamf_inds = (dd-1)*length( chamfer_k )+1:dd*length( chamfer_k )
            chamfer_dists[ rep, : ] = c_ii
        end
        if eCDF == 1
            # cdfs_inds = (dd-1)*data[:nbin]*length( LL )+1:dd*data[:nbin]*length( LL )
            cdfs[ rep, : ] = cdf_ii
        end

        # Calculate summaries for muu and C based on simulated data
        if case == "bsl"
            cdf_ii, c_ii = create_summaries( x, yhat, data, create_kdtree, kn )
            if chamfer == 1
                # chamf_inds = (dd-1)*length( chamfer_k )+1:dd*length( chamfer_k )
                chamfer_dists_hat[ rep, : ] = c_ii
            end
            if eCDF == 1
                # cdfs_inds = (dd-1)*data[:nbin]*length( LL )+1:dd*data[:nbin]*length( LL )
                cdfs_hat[ rep, : ] = cdf_ii
            end

            data[:binss] = Vector{Vector{Float64}}(undef, nL)
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

        if any(isnan, C)
            println(summary_stats_hat)
            sleep(3)
        end
    end

    return data, summary_stats
end

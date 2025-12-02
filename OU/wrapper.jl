# For cases where difference and no 2D data is used
function Wrapper(theta::Vector{Float64}, data::Dict{Symbol, Any})
    dt = data[:synth_dt]
    Ndata = data[:Ndata]
    init = data[:synth_init]
    R0_all = data[:R0_all]
    nobs = data[:nobs]
    case_nsim, case_nrep = data[:case_dim]
    chamfer = data[:chamfer]
    chamfer_k = data[:chamfer_k]
    eCDF = data[:eCDF]
    nL = data[:nL]
    diff_order = data[:diff_order]
    use_diff = data[:use_diff]
    bins = data[:bins]
    LL = data[:LL]


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
        X = OU_solve( theta, init, Ndata; dt = dt )
        R_sim_all[ sim, : ] = X
        if use_diff == 1
            X_diffs = calculate_diffs( X, diff_order, dt )
            for ii in eachindex( diff_order )
                R_sim_diff_all[ii][ sim, : ] = X_diffs[ii]
            end
        end
    end

    if use_diff == 1
        R_sim_all = [ R_sim_all, R_sim_diff_all... ]  # Concatenate all simulated data into vector of matrices
    else
        R_sim_all = [ R_sim_all ]
    end

    cdfs_ii_2 = Vector{Float64}(undef, 0)
    chamfer_ii = Vector{Float64}(undef, 0)
    # Create eCDFs and chamfer distances for each simulations
    for dd in eachindex(R0_all) # Loop over data types: original and differences
        x = R0_all[dd]
        R_sim_dd = R_sim_all[dd]
        data[:bins] = bins[dd]

        for (ii, y) in enumerate( eachrow( R_sim_dd ) )
            if -1 in data[:LL]          # Signal feature
                cdfs_ii, _ = empcdf( y, nx=data[:nbin], x=data[:bins][1] )
                data[:LL] = data[:LL][2:end]
                data[:bins] = data[:bins][2:end]
                if ~isempty( data[:LL] ) || chamfer == 1
                    if data[:resample] != 0
                        data, cdfs_ii_2, chamfer_ii = resample_data( x, y, data )
                    end

                end
                data[:LL] = LL
                data[:bins] = bins[dd]

                cdfs_ii = vcat( vec(cdfs_ii), cdfs_ii_2 )
            else
                if data[:resample] != 0
                    cdfs_ii, chamfer_ii = resample_data( y, data )
                end
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

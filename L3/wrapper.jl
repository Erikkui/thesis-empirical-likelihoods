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
    diff_order = data[:diff_order]
    use_diff = data[:use_diff]
    bins = data[:bins]


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
    R_sim_all = Vector{ Matrix{Float64} }(undef, case_nsim)
    if use_diff == 1
        R_sim_diff_all = Vector{ Vector{ Matrix{Float64} } }(undef, length( diff_order ) )
        for ii in eachindex( diff_order )
            R_sim_diff_all[ii] = Vector{ Matrix{Float64} }(undef, case_nsim)
        end
    end
    for sim in 1:case_nsim
        N = lorenz_solve( init, theta, N_end; dt = dt )
        R_sim_all[ sim ] = N

        if use_diff == 1
            N_diffs = calculate_diffs( N, diff_order, dt )
            for ii in eachindex( diff_order )
                R_sim_diff_all[ii][sim] = @view N_diffs[ii][ :, : ]
            end
        end
    end


    if use_diff == 1
        R_sim_all = [ R_sim_all, R_sim_diff_all... ]  # Concatenate all simulated data into vector of vector of matrices
    else

        R_sim_all = [ R_sim_all ]
    end

    # For each data set: signal and its difference(s)
    for dd in eachindex(R0_all)
        x = R0_all[dd]
        Rsim = R_sim_all[dd]
        data[:bins] = bins[dd]

        for (ii, y) in enumerate( Rsim )

            cdfs_ii, chamfer_ii = create_summaries( x, y, data, create_kdtree, kn )

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

    println(summary_stats)

    if data[:log] == "log"
        summary_stats = log.(summary_stats)
    end
    # sleep(1)

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

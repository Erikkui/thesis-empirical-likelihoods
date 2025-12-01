function resample_data( R_all, data::Dict{Symbol, Any} )
    # Create likelihood by resampling data in R0. Create bins if using ecdfs.
    # R0: nepo sets of nobs vectors each, pooled here together for resampling
    # Either eCDF, Chamfer distance, or both returned in 'cdf'
    # For eCDF the bin values are returned in 'data', computed by the
    # 'nsamp' values of the first inner loop.
    #
    # data.res_dim =[method,nrepo,nsamp]:
    # nrep      size of outer loop
    # nsamp     size of inner loop

    nrep, nsamp = data[:res_dim]
    eCDF = data[:eCDF]
    chamfer = data[:chamfer]

    create_kdtree = data[:create_kdtree]
    kn = data[:kn]

    ntot = size( R_all, 2 )  # Total number of observations

    # Initialize data structures
    chamfer_dists = Array{Float64}(undef, nrep*nsamp, 0)
    cdfs = Array{Float64}(undef, nrep*nsamp, 0)
    bins = nothing
    if chamfer == 1
        chamfer_k = data[:chamfer_k]
        chamfer_dists = zeros( nrep*nsamp, length( chamfer_k ) )
    end
    if eCDF == 1
        nL = data[:nL]
        LL = data[:LL]
        nbin = data[:nbin]

        cdfs = zeros( nrep*nsamp, nL*nbin )

        # For bin limits and calculation
        minmax = [ fill(-Inf, (1, nL) ); fill(Inf, (1, nL)) ]
        mins = copy( minmax[1, :] )
        maxs = copy( minmax[2, :] )
    end

    # Setting iterator. Use progress bar if calculating summaries in addition to bins
    if nrep > 1
        iter = ProgressBar( 1:nrep*nsamp )
    else
        iter = 1:nrep*nsamp
    end

    for ii in iter
        # OOB sampling
        in_bag = rand( 1:ntot, ntot)
        unique!( in_bag )
        out_of_bag = setdiff( 1:ntot, in_bag )
        x = @view R_all[ :, in_bag ]
        y = @view R_all[ :, out_of_bag ]

        # Find bin limits and calculate bins OR only calculate summaries for this sample
        bins = Vector{ Vector{Float64} }( undef, 0 )
        if ii <= nsamp && eCDF == 1

            cdf_ii, c_ii = create_summaries( x, y, data, create_kdtree, kn; same = true )

            # Min and max values for bin creation
            for jj in 1:nL
                m1 = minimum( cdf_ii[jj] )
                mins[jj] = max( mins[jj], m1 )
                M1 = maximum( cdf_ii[jj] )
                maxs[jj] = min( maxs[jj], M1 )
            end

            if ii == nsamp  # Only at the nsamp'th sample, compute bins

                minmax[1, :] = mins
                minmax[2, :] = maxs
                # data[:minmax] = minmax

                for jj in 1:nL
                    data[:minmax] = minmax[:, jj]
                    if data[:uni] == "yax"
                        data[:r] = cdf_ii[jj]    # data[:r][LL] = signal/CIL/ID feature vector
                    end
                    bins_nL = bin_select( data )
                    push!( bins, bins_nL )
                end
            end

        else    # Only calculate summaries for this sample
            cdf_ii, c_ii = create_summaries( x, y, data, create_kdtree, kn )
            cdfs[ ii, : ] = cdf_ii
            chamfer_dists[ ii, : ] = c_ii
        end

        # For progress bar description
        if nrep > 1
            set_description(iter, "Resampling data:")
        end
    end

    summary_stats = [ cdfs chamfer_dists ]

    # Remove rows of all zeros
    mask = any( abs.( summary_stats ) .> 10e-16, dims=2 )
    summary_stats = summary_stats[ vec(mask), : ]

    if data[:log] == "log"
        summary_stats = log.(summary_stats)
    end

    return data, summary_stats, bins
end

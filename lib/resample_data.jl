function resample_bins( R_all, data::Dict{Symbol, Any} )
    # Resample data to create bins for eCDF calculation
    nsamp = data[:res_dim][2]
    create_kdtree = data[:create_kdtree]
    kn = data[:kn]
    nL = data[:nL]

    # For bin limits and calculation
    minmax = zeros( 2, nL )
    mins = zeros( nsamp, nL)
    maxs = zeros( nsamp, nL)

    ntot = size( R_all, 2 )  # Total number of observations
    bins = Vector{ Vector{Float64} }( undef, 0 )
    for ii in 1:nsamp
        # OOB sampling
        in_bag = rand( 1:ntot, ntot)
        unique!( in_bag )
        out_of_bag = setdiff( 1:ntot, in_bag )
        x = copy( R_all[ :, in_bag ] )
        y = copy( R_all[ :, out_of_bag ] )

        cdf_ii, _ = create_summaries( x, y, data, create_kdtree, kn; same = true )

        # Min and max values for bin creation
        m = minimum.( cdf_ii )
        M = maximum.( cdf_ii )
        mins[ii, :] = m
        maxs[ii, :] = M

        if ii == nsamp  # Only at the nsamp'th sample, compute bins

            minmax[1, :] = maximum( mins, dims=1 )
            minmax[2, :] = minimum( maxs, dims=1 )
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
    end
    return data, bins
end


function resample_data( R_all, data::Dict{Symbol, Any} )
    ntot = size( R_all, 2 )  # Total number of observations
    nrep = data[:res_dim][1]
    chamfer = data[:chamfer]
    eCDF = data[:eCDF]
    create_kdtree = data[:create_kdtree]
    kn = data[:kn]

    # Resample data to create data covariance/mean
    chamfer_dists = Array{Float64}(undef, nrep, 0)
    cdfs = Array{Float64}(undef, nrep, 0)

    if chamfer == 1
        chamfer_k = data[:chamfer_k]
        chamfer_dists = zeros( nrep, length( chamfer_k ) )
    end
    if eCDF == 1
        nL = length( data[:LL] )
        nbin = data[:nbin]
        cdfs = zeros( nrep, nL*nbin )
    end

    for ii in 1:nrep
        # OOB sampling
        in_bag = rand( 1:ntot, ntot)
        unique!( in_bag )
        out_of_bag = setdiff( 1:ntot, in_bag )
        x = copy( R_all[ :, in_bag ] )
        y = copy( R_all[ :, out_of_bag ] )

        cdf_ii, c_ii = create_summaries( x, y, data, create_kdtree, kn )
        cdfs[ ii, : ] = cdf_ii
        chamfer_dists[ ii, : ] = c_ii
    end

    return data, cdfs, chamfer_dists
end


function resample_data( R_obs, R_sim, data::Dict{Symbol, Any} )
    nrep = data[:case_dim][2]
    chamfer = data[:chamfer]
    eCDF = data[:eCDF]
    create_kdtree = data[:create_kdtree]
    kn = data[:kn]
    ntot = size( R_sim, 2 )  # Total number of observations

    # Resample data to create data covariance/mean
    chamfer_dists = Array{Float64}(undef, nrep, 0)
    ecdfs = Array{Float64}(undef, nrep, 0)

    if chamfer == 1
        chamfer_k = data[:chamfer_k]
        chamfer_dists = zeros( nrep, length( chamfer_k ) )
    end
    if eCDF == 1
        nL = length( data[:LL] )
        nbin = data[:nbin]
        ecdfs = zeros( nrep, nL*nbin )
    end

    for ii in 1:nrep
        # OOB sampling
        in_bag = rand( 1:ntot, ntot)
        unique!( in_bag )
        y = copy( R_sim[:, in_bag] )

        cdf_ii, c_ii = create_summaries( R_obs, y, data, create_kdtree, kn )
        ecdfs[ ii, : ] = cdf_ii
        chamfer_dists[ ii, : ] = c_ii
    end

    # if size( ecdfs, 1 ) > 1
    #     ecdfs = mean( ecdfs, dims=1 )
    #     chamfer_dists = mean( chamfer_dists, dims=1 )
    # end

    return data, ecdfs, chamfer_dists
end

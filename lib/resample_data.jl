function resample_bins( R_all, data::Dict{Symbol, Any} )
    # Resample data to create bins for eCDF calculation
    res_bins = data[:res_dim][2]
    create_kdtree = data[:create_kdtree]
    kn = data[:kn]
    nL = data[:nL]

    # For bin limits and calculation
    minmax = zeros( 2, nL )
    mins = zeros( res_bins, nL)
    maxs = zeros( res_bins, nL)

    ntot = size( R_all, 2 )  # Total number of observations
    bins = Vector{ Vector{Float64} }( undef, 0 )
    for ii in 1:res_bins

        # If signal, resample a continuous partition. Else OOB sampling
        if -1 in data[:LL]
            window = data[:window]
            start_ind = rand( 1:ntot-window )
            end_ind = start_ind+window
            y = copy( R_all[ 1:1, start_ind:end_ind] )
            x = copy( R_all[ 1:1, 1:kn+1 ] )    # Not actually needed in bin creation
        else
            # OOB sampling
            in_bag = rand( 1:ntot, ntot)
            unique!( in_bag )
            out_of_bag = setdiff( 1:ntot, in_bag )
            x = copy( R_all[ :, in_bag ] )
            y = copy( R_all[ :, out_of_bag ] )
        end

        cdf_ii, _ = create_summaries( x, y, data, create_kdtree, kn )

        # Min and max values for bin creation
        m = minimum.( cdf_ii )
        M = maximum.( cdf_ii )
        mins[ii, :] = m
        maxs[ii, :] = M

        if ii == res_bins  # Only at the nsamp'th sample, compute bins

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

    # Resample data to create data mean
    chamfer_dists = Array{Float64}(undef, nrep, 0)
    cdfs = Array{Float64}(undef, nrep, 0)

    if chamfer == 1
        chamfer_k = data[:chamfer_k]
        chamfer_dists = zeros( nrep, length( chamfer_k ) )
    end
    if eCDF == 1
        nL = data[:nL]
        nbin = data[:nbin]
        cdfs = zeros( nrep, nL*nbin )
    end

    for ii in 1:nrep
        # If signal, resample a continuous partition. Else OOB sampling
        if -1 in data[:LL]
            window = data[:window]
            start_ind = rand( 1:length(R_all)-window )
            y = copy( R_all[ 1:1, start_ind:start_ind+window] )
            x = copy( R_all[ 1:1, 1:kn+1 ] )    # Not actually needed in data mean calculation
        else
            # OOB sampling
            in_bag = rand( 1:ntot, ntot)
            unique!( in_bag )
            out_of_bag = setdiff( 1:ntot, in_bag )
            x = copy( R_all[ :, in_bag ] )
            y = copy( R_all[ :, out_of_bag ] )
        end

        cdf_ii, c_ii = create_summaries( x, y, data, create_kdtree, kn )
        cdfs[ ii, : ] = cdf_ii
        chamfer_dists[ ii, : ] = c_ii
    end

    return data, cdfs, chamfer_dists
end


function resample_data( R_data, R_sim, data::Dict{Symbol, Any} )
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
        if -1 in data[:LL]
            cdf_ii, c_ii = create_summaries( R_data, R_sim, data, create_kdtree, kn )
        else
            # OOB sampling
            in_bag = rand( 1:ntot, ntot)
            unique!( in_bag )
            y = copy( R_sim[:, in_bag] )
            cdf_ii, c_ii = create_summaries( R_data, y, data, create_kdtree, kn )
        end
        ecdfs[ ii, : ] = cdf_ii
        chamfer_dists[ ii, : ] = c_ii
    end

    if size( ecdfs, 1 ) > 1
        ecdfs = mean( ecdfs, dims=1 )
        chamfer_dists = mean( chamfer_dists, dims=1 )
    end

    return data, ecdfs, chamfer_dists
end

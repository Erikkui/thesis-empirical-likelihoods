# When called in GSL case for the first time
function create_likelihood( data::Dict{Symbol, Any} )
    # Create likelihood either by resampling or by epochs

    # Check if kdtree needs to be created: either chamfer = 1 or ID features are used
    data = determine_kdtree( data )

    # Concatenate all nepo R0 matrices, R0 = Vector{ Matrix{Float64} }
    data[:R0_all] = [ hcat( data[:R0]... ) ]    # Vector a single dim x nobs*nepo matrix

    # Initialize bin storage
    data[:bins] = Vector{ Any }( nothing, length( data[:R0_all] ) )
    if data[:eCDF] == 1
        data[:bins] = Vector{ Vector{ Vector{Float64} } }(undef, 1 + length(data[:diff_order])*data[:use_diff] )
    end

    if data[:use_diff] == 1
        data[:R0_diff_all] = [ hcat( data[:R0_diff][ii]... ) for ii in eachindex(data[:R0_diff]) ]
        data[:R0_all] = [ data[:R0_all]..., data[:R0_diff_all]... ]
        # Vector of all matrices: first is R0, then R0_diff for each diff order.
        # Size of each matrix is dim x (nobs*nepo) (note that nobs is different for each diff order)
    end

    data[:bins_done] = false
    # Create bins: bins[ii] is for signal/difference, and each each row is for CIL/ID feature
    if data[:eCDF] == 1
        for ii in eachindex( data[:R0_all] )
            R0_ii = data[:R0_all][ii]
            data, bins = resample_bins( R0_ii, data )
            data[:bins][ii] = bins
        end
    end
    data[:bins_done] = true

    LL = deepcopy( data[:LL] )  # Store original LL to restore later
    bins = deepcopy( data[:bins] )  # Store bins to restore later
    ecdfs = Vector{Matrix{Float64}}(undef, 0)
    chamfer_dists = Vector{Matrix{Float64}}(undef, 0)
    for ii in eachindex( data[:R0_all] )    # For each signal/difference
        RR_ii = data[:R0_all][ii]
        data[:bins] = bins[ii]
        data, cdfs_ii, chamfer_ii = resample_data( RR_ii, data )

        # Normalize chamfer due to being possibly very large compared to cdf values
        means = mean(chamfer_ii, dims = 1)
        stds = std( chamfer_ii, dims = 1)
        chamfer_ii = (chamfer_ii .- means) ./ stds

        push!( ecdfs, cdfs_ii )
        push!( chamfer_dists, chamfer_ii )
    end
    ecdfs = hcat( ecdfs... )
    chamfer_dists = hcat( chamfer_dists... )
    summary_stats = hcat( ecdfs, chamfer_dists )

    data[:bins] = bins
    data[:LL] = LL

    if data[:log] == "log"
        summary_stats = log.( summary_stats )
    end

    # Compute mean
    muu = mean( summary_stats, dims=1 )

    # Compute covariance
    if data[:C_how] == "don"
        C = donsker( muu, data[:nobs] )
    elseif data[:C_how] == "cov"
        C = cov( summary_stats )
    end

    data[:muu_data] = muu
    data[:C] = C

    return data, summary_stats
end


function determine_kdtree( data )
    kn_c = 0
    kn_e = 0
    if data[:chamfer] == 1
        kn_c = maximum( data[:chamfer_k] )
    end
    if data[:eCDF] == 1
        kn_e = maximum( data[:LL] )
    end
    kn = max( kn_c, kn_e )+1  # Number of neighbors to find if needed
    kdtree_condition = ( kn_e > 0 || data[:chamfer] == 1 )
    create_kdtree = kdtree_condition ? true : false
    data[:create_kdtree] = create_kdtree
    data[:kn] = kn

    return data
end

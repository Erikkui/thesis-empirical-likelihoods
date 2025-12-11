
function create_summaries( x, y, data, make_kdtree, kn )
    # Create summary statistics: ecdf's and/or chamfer distances

    # Whether to build KDTree
    D_xy =  Vector{Vector{Float64}}(undef, 0)
    cdfs = Vector{Vector{Float64}}(undef, 0)
    chamfers = Vector{Float64}(undef, 0)
    if make_kdtree
        ytree = KDTree(y)
        _, D_xy = knn( ytree, x, kn, true )
        data[:D_xy] = D_xy
    end

    # Calculate full pairwise matrix (vector):
    if 0 in data[:LL]
        D = pairwise_turbo( x, y )
    end

    # If ecdf's are used
    if data[:eCDF] == 1
        if -1 in data[:LL]          # Signal feature
            cdfs = cdf_do( y, data )
        else
            cdfs = cdf_do( D, data )
        end
    end

    # If chamfer distances are used
    if data[:chamfer] == 1
        if make_kdtree
            chamfers = chamfer_dist( x, y, D_xy, data[:chamfer_k] )
        else
            chamfers = chamfer_dist( D, data[:D_xy_sorted], data[:chamfer_k] )
        end
    end

    # If bins have been calculated, return cdfs as one long vector; otherwise return them
    # as vector of vectors from which bins can be calculated for each feature.
    if data[:bins_done]
        return vcat( cdfs... ), chamfers
    else
        return cdfs, chamfers
    end
end


# Special case for -1
function create_summaries( D, data )
    # Create summary statistics: ecdf's from signals only
    c_ii = Vector{Float64}(undef, 0)
    cdf_ii = Vector{Vector{Float64}}(undef, 0)
    cdf_ii_temp = cdf_do( D, data )
    push!( cdf_ii, cdf_ii_temp... )
    # If bins have been calculated, return cdfs as one long vector; otherwise return them
    # as vector of vectors from which bins can be calculated for each feature.
    if data[:bins_done]
        return vcat( cdf_ii... ), c_ii
    else
        return cdf_ii, c_ii
    end
end

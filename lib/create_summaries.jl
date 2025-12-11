
function create_summaries( x, y, data, make_kdtree, kn; same = false )
    # Create summary statistics: ecdf's and/or chamfer distances

    # Whether to build KDTree
    D_xy =  Vector{Vector{Float64}}(undef, 0)
    cdf_ii = Vector{Vector{Float64}}(undef, 0)
    c_ii = Vector{Float64}(undef, 0)
    if make_kdtree
        ytree = KDTree(y)
        _, D_xy = knn( ytree, x, kn, true )
        data[:D_xy] = D_xy
    end

    if 0 in data[:LL] || ( !make_kdtree && ( data[:chamfer] == 1 || any( data[:LL] .>= 0 ) ) )
        D = pairwise_turbo( x, y )

        if ( data[:chamfer] == 1 || any( data[:LL] .>= 1 ) ) && !make_kdtree
            D_xy_sorted = Matrix{Float64}( undef, size(x, 2), kn )
            for ii in axes(D, 1)
                row_view = view(D, ii, :)
                D_xy_sorted[ii, :] = partialsort( row_view, 1:kn )
            end
            data[:D_xy_sorted] = D_xy_sorted
        end
    end


    # If ecdf's are used
    if data[:eCDF] == 1
        if -1 in data[:LL]          # Signal feature
            if same
                cdf_ii_temp = cdf_do( x, data )
            else
                XY = hcat( x, y )
                cdf_ii_temp = cdf_do( XY, data )
            end
            push!( cdf_ii, cdf_ii_temp... )
        end
        if any( data[:LL] .>= 0 )   # Distance (CIL/ID) features
            cdf_ii_temp = cdf_do( D, data )
            push!( cdf_ii, cdf_ii_temp... )
        end

    end

    # If chamfer distances are used
    if data[:chamfer] == 1
        if make_kdtree
            c_ii = chamfer_dist( x, y, D_xy, data[:chamfer_k] )
        else
            c_ii = chamfer_dist( D, data[:D_xy_sorted], data[:chamfer_k] )
        end
    end

    # If bins have been calculated, return cdfs as one long vector; otherwise return them
    # as vector of vectors from which bins can be calculated for each feature.
    if data[:bins_done]
        return vcat( cdf_ii... ), c_ii
    else
        return cdf_ii, c_ii
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


function create_summaries( x, y, data, make_kdtree, kn )
    # Create summary statistics: ecdf's and/or chamfer distances

    # Whether to build KDTree: always if chamfer used, or if ID features are used
    D_xy =  Vector{Vector{Float64}}(undef, 0)
    cdf_ii = Vector{Vector{Float64}}(undef, 0)
    c_ii = Vector{Float64}(undef, 0)
    if make_kdtree
        ytree = KDTree(y)
        _, D_xy = knn( ytree, x, kn, true )
    end

    # If ecdf's are used
    nL = data[:nL]
    if data[:eCDF] == 1
        D = Matrix{Float64}(undef, 0, 0)
        if 0 in data[:LL]
            if -2 in data[:LL]      # If -2, then do not use distances
                D = hcat( x, y )
                cdf_ii_temp = cdf_do( D, data )
                push!( cdf_ii, cdf_ii_temp... )

            elseif -1 in data[:LL]  # If -1, then use distances and signal
                bins_all = copy( data[:binss] )
                data[:binss] = bins_all[1:nL-1]  # Pick bins for signal

                D = hcat( x, y )
                cdf_ii_temp = cdf_do( D, data ) # cdf for signal
                push!( cdf_ii, cdf_ii_temp... )

                data[:D_xy] = D_xy
                D = pairwise( Euclidean(), x, y )
                data[:binss] = bins_all[nL:end]  # Pick bins for distances

                cdf_ii_temp = cdf_do( D, data ) # cdf for distances
                push!( cdf_ii, cdf_ii_temp... )

                data[:binss] = bins_all     # Restore bins to data
                pop!( data, :D_xy )  # Remove D_xy from data
            else                    # Else: only use distances
                D = pairwise( Euclidean(), x, y )
                cdf_ii_temp = cdf_do( D, data )
                push!( cdf_ii, cdf_ii_temp... )
            end
        else

            D = pairwise( Euclidean(), x, y )
            cdf_ii_temp = cdf_do( D, data )
            push!( cdf_ii, cdf_ii_temp... )
        end
    end

    # If chamfer distances are used
    if data[:chamfer] == 1
        c_ii = chamfer_dist( x, y, D_xy, data[:chamfer_k] )
    end

    # If bins have been calculated, return cdfs as one long vector; otherwise return them
    # as vector of vectors from which bins can be calculated for each feature.
    if isassigned( data[:binss], 1 )
        return vcat( cdf_ii... ), c_ii
    else
        return cdf_ii, c_ii
    end
end


# Special case for -2
function create_summaries( D, data )
    # Create summary statistics: ecdf's from signals only
    c_ii = Vector{Float64}(undef, 0)
    cdf_ii = Vector{Vector{Float64}}(undef, 0)
    cdf_ii_temp = cdf_do( D, data )
    push!( cdf_ii, cdf_ii_temp... )

    # If bins have been calculated, return cdfs as one long vector; otherwise return them
    # as vector of vectors from which bins can be calculated for each feature.
    if isassigned( data[:binss], 1 )
        return vcat( cdf_ii... ), c_ii
    else
        return cdf_ii, c_ii
    end
end

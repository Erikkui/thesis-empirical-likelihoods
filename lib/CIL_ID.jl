# Standard use case
function CIL_ID( distances, LL::Vector{Int} )
    # Calculate CIL and/or ID vectors from the pairwise distances

    N_features = sum(LL .>= 0 )
    features = Vector{Vector{Float64}}( undef, N_features )
    feature_ind = 1

    if 0 in LL  # CIL feature
        cil_dists = distances[:]
        features[ feature_ind ] = cil_dists
        filter!( >(0), LL )
        feature_ind += 1
    end

    if ~isempty(LL)     # ID features
        n_dist_rows = size( distances, 1 )
        n_ID_features = length(LL)
        ID_max = maximum(LL)
        ID_ratios = zeros( Float64, n_dist_rows, ID_max )

        for ii in 1:n_dist_rows     # ID ratios
            dists_row = @view distances[ ii, : ]
            println(dists_row)
            partialsort!( dists_row, 1:ID_max+1 )   # Sort the up to ID_max+1 elements

            ratios_ii = @views( dists_row[ 2:ID_max+1 ] ./ dists_row[ 1:ID_max ] )
            ID_ratios[ ii, : ] = ratios_ii
        end

        for ii in 1:n_ID_features
            ID_ind = LL[ii]
            features[ feature_ind ] = @view ID_ratios[ :, ID_ind ]
            feature_ind += 1
        end

    end

    return features
end



# For kdtree cases
function CIL_ID( D, D_xy::Vector{Vector{Float64}}, LL::Vector{Int} )

    N_features = sum(LL .>= 0 )
    features = Vector{Vector{Float64}}( undef, N_features )
    feature_ind = 1

    if 0 in LL  # CIL feature
        cil_dists = D[:]
        features[1] = cil_dists
        feature_ind = 2
    end
    filter!( >(0), LL )

    # ID feature
    if ~isempty(LL)
        n_dist_rows = length( D_xy )
        n_ID_features = length(LL)
        ID_max = maximum(LL)
        ID_ratios = zeros( Float64, n_dist_rows, ID_max )

        for ii in 1:n_dist_rows     # ID ratios
            dists_row = D_xy[ii]
            ratios_ii = dists_row[ 2:ID_max+1 ] ./ dists_row[ 1:ID_max ]

            ID_ratios[ ii, : ] = ratios_ii
        end

        for ii in 1:n_ID_features
            ID_ind = LL[ii]
            features[ feature_ind ] = ID_ratios[ :, ID_ind ]
            feature_ind += 1
        end
    end

    return features

end

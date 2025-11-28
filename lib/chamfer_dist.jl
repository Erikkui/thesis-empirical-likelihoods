# Non-kdtree version
function chamfer_dist( D, chamfer_dims::Vector{Int} )

    n, m = size(D)

    if length( chamfer_dims ) == 1 && chamfer_dims[1] == 1  # Standard chamfer distance
        x_to_y = minimum( D, dims=2 )
        y_to_x = minimum( D, dims=1 )
        chamfer_dist = sum( x_to_y )/n + sum( y_to_x )/m
    else
        # Sort the distance matrix
        Dsorted_xy = sort( D, dims=2 )
        Dsorted_yx = sort( D, dims=1 )

        # Extract the chamfer dimensions
        x_to_y = @view( Dsorted_xy[ :, chamfer_dims ] )
        y_to_x = @view( Dsorted_yx[ chamfer_dims, : ] )

        # Compute the chamfer distance
        chamfer_dist = sum( x_to_y, dims = 1 )./n + sum( y_to_x, dims = 2 )'./m

    end

    return chamfer_dist
end


# kdtree version
function chamfer_dist( x, y, D_xy::Vector{Vector{Float64}}, chamfer_dims::Vector{Int} )
    # Compute Chamfer distance between x and y
    n = size( x, 2 )
    m = size( y, 2 )
    k = maximum( chamfer_dims )

    tree_x = KDTree( x )
    _, D_yx = knn( tree_x, y, k, true )

    # Extract the chamfer dimensions
    x_to_y = @view sum( D_xy )[ chamfer_dims ]
    y_to_x = @view sum( D_yx )[ chamfer_dims ]

    # Compute the chamfer distance
    chamfer_dist = x_to_y./n + y_to_x./m

    return chamfer_dist
end

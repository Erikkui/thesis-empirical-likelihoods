function chamfer_dist(D, neighbors)
    kmax = maximum( neighbors )
    N_xy, N_yx = size(D)

    # Nearest neighbors for each x
    dists_x_to_y = sum( @views( D[:, 1:kmax] ) , dims=1 )

    # Nearest neighbors for each y
    sort!( D; dims=1 )
    dists_y_to_x = sum( @views( D[1:kmax, :] ) , dims=2 )'

    # Compute chamfer distance
    chamfer_dist = dists_x_to_y./N_xy + dists_y_to_x./N_yx

    return chamfer_dist
end


# kdtree version
function chamfer_dist( x, y, D_xy::Vector{Vector{Float64}}, chamfer_dims::Vector{Int} )
    # Compute Chamfer distance between x and y
    n = size( x, 2 )
    m = size( y, 2 )
    k = maximum( chamfer_dims )

    tree_x = KDTree( x )

    # println( size(x), "  ", size(y), "\n\n", k )
    _, D_yx = knn( tree_x, y, k, false )

    # Extract the chamfer dimensions
    x_to_y = @view sum( D_xy )[ chamfer_dims ]
    y_to_x = @view sum( D_yx )[ chamfer_dims ]

    # Compute the chamfer distance
    chamfer_dist = x_to_y./n + y_to_x./m

    return chamfer_dist
end

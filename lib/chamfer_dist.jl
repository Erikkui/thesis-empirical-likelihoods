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

    _, D_yx = knn( tree_x, y, k, true )

    # Extract the chamfer dimensions
    x_to_y = @view sum( D_xy )[ chamfer_dims ]
    y_to_x = @view sum( D_yx )[ chamfer_dims ]

    # Compute the chamfer distance
    chamfer_dist = x_to_y./n + y_to_x./m

    return chamfer_dist
end

function chamfer_dist(D, D_sorted, neighbors)
    kmax = maximum( neighbors )
    N_xy, N_yx = size(D)

    # Nearest neighbors for each x
    dists_x_to_y = sum( @views( D_sorted[:, 1:kmax] ) , dims=1 )

    # Nearest neighbors for each y
    D_sorted = Matrix{Float64}( undef, kmax, N_yx )
    for ii in axes(D, 2)
        col_view = view(D, :, ii)
        D_sorted[:, ii] = partialsort( col_view, 1:kmax )
    end
    dists_y_to_x =  sum( D_sorted , dims=2 )'

    # Compute chamfer distance
    chamfer_dist = dists_x_to_y./N_xy + dists_y_to_x./N_yx

    return chamfer_dist
end

function cdf_do( D, data::Dict{Symbol, Any} )
    LL = copy( data[:LL] )
    nL = data[:nL]

    # If CIL/ID is calculated
    if any( LL .>= 0 )
        D_xy = data[:D_xy]
        features = CIL_ID( D, D_xy, LL )
    else
        features = [ vec(D) ]
    end

    # If bins are found, calculate ecdfs and return them; otherwise only return features
    if data[:bins_done]
        bins = data[:bins]
        nbin = data[:nbin]

        # Compute CDFs for given bins
        cdfs = Vector{Vector{Float64}}(undef, 0)

        for L in 1:nL
            cdf_L, _ = empcdf( features[L], nx=nbin, x=bins[L] )
            push!(cdfs, vec(cdf_L))
        end

        return cdfs
    else
        return features
    end

end

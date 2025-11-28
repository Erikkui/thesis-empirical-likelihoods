# When called in GSL case for the first time
function create_likelihood( data::Dict{Symbol, Any} )
    # Create likelihood either by resampling or by epochs

    # Check if signal-only likelihood is to be used
    if -2 in data[:LL]
        data[:signal_only] = true
    end

    # Check if kdtree needs to be created: either chamfer = 1 or ID features are used
    kn_c = 0
    kn_e = 0
    if data[:chamfer] == 1
        kn_c = length( data[:chamfer_k] )
    end
    if data[:eCDF] == 1
        kn_e = maximum( data[:LL] )
    end
    kn = max( kn_c, kn_e )  # Number of neighbors to find if needed
    create_kdtree = kn_e > 0 || data[:chamfer] == 1 ? true : false
    data[:create_kdtree] = create_kdtree
    data[:kn] = kn

    if data[:resample] == 1     # Make bins and CDF/chamfer dist by resampling
        # If using resampling, store all R0 in one big matrix

        # Concatenate all nepo R0 matrices
        data[:R0_all] = [ hcat( data[:R0]... ) ]

        if data[:use_diff] == 1
            data[:R0_diff_all] = [ hcat( data[:R0_diff][ii]... ) for ii in eachindex(data[:R0_diff]) ]
            data[:R0_all] = [ data[:R0_all]..., data[:R0_diff_all]... ]
        end

        if data[:case] == "gsl"
            bins_all = Vector{Vector{Float64}}(undef, 0)
            summary_stats = Vector{ Matrix{Float64} }(undef, 0)

            for ii in eachindex( data[:R0_all] )
                data[:binss] = Vector{Vector{Float64}}(undef, 1)
                data, summary_stats_ii, _ = resample_data( data[:R0_all][ii], data )
                push!( summary_stats, summary_stats_ii )
                push!( bins_all, data[:binss]... )
            end
            summary_stats = hcat( summary_stats... )
            data[:binss] = bins_all
        else
            # If using BSL, likelihood is computed according to simulations
            return data, nothing
        end

    else
        # TODO
        # Make bins by epochs
        # data, cdf = pairwise_distL(data[:R0], data)
    end

    # Compute mean
    if size(summary_stats, 1) > 1
        muu = mean( summary_stats, dims=1 )
    else
        muu = summary_stats
    end

    # Compute covariance
    if data[:C_how] == "don"
        C = donsker( muu, data[:nobs] )
    elseif data[:C_how] == "cov"
        C = cov( summary_stats )
    end

    data[:muu] = muu
    data[:C] = C

    return data, summary_stats
end

# When called in GSL case for the first time
function create_likelihood( data::Dict{Symbol, Any} )
    # Create likelihood either by resampling or by epochs

    # Check if kdtree needs to be created: either chamfer = 1 or ID features are used
    kn_c = 0
    kn_e = 0
    if data[:chamfer] == 1
        kn_c = length( data[:chamfer_k] )
    end
    if data[:eCDF] == 1
        kn_e = maximum( data[:LL] )
    end
    kn = max( kn_c, kn_e )+1  # Number of neighbors to find if needed
    kdtree_condition = ( kn_e > 0 || data[:chamfer] == 1 ) && data[:data_dim] > 1
    create_kdtree = kdtree_condition ? true : false
    data[:create_kdtree] = create_kdtree
    data[:kn] = kn


    # Initialize bin storage
    data[:bins] = Dict{Symbol, Any}()
    if data[:eCDF] == 1
        data[:bins] = Vector{ Vector{ Vector{Float64} } }(undef, 1 + length(data[:diff_order])*data[:use_diff] )
    end

    if data[:resample] != 0     # Make bins and CDF/chamfer dist by resampling
        # If using resampling, store all R0 in one big matrix

        # Concatenate all nepo R0 matrices, R0 = Vector{ Matrix{Float64} }
        data[:R0_all] = [ hcat( data[:R0]... ) ]    # Vector a single dim x nobs*nepo matrix

        if data[:use_diff] == 1
            data[:R0_diff_all] = [ hcat( data[:R0_diff][ii]... ) for ii in eachindex(data[:R0_diff]) ]
            data[:R0_all] = [ data[:R0_all]..., data[:R0_diff_all]... ]
            # Vector of all matrices: first is R0, then R0_diff for each diff order.
            # Size of each matrix is dim x (nobs*nepo) (note that nobs is different for each diff order)
        end

        data[:bins_done] = false
        if data[:case] == "gsl"
            #TODO maybe to be deprecated
        else
            # Create bins
            for ii in eachindex( data[:R0_all] )
                R0_ii = data[:R0_all][ii]
                data, bins = resample_bins( R0_ii, data )
                data[:bins][ii] = bins
            end
            data[:bins_done] = true

            LL = data[:LL]  # Store original LL
            bins = data[:bins]  # Store bins for later use
            ecdfs = Vector{Matrix{Float64}}(undef, 0)
            chamfer_dists = Vector{Matrix{Float64}}(undef, 0)
            for ii in eachindex( data[:R0_all] )
                RR_ii = data[:R0_all][ii]
                data[:bins] = bins[ii]

                if -1 in data[:LL]          # Signal feature
                    cdfs_ii, _ = empcdf( vec(RR_ii), nx=data[:nbin], x=data[:bins][1] )
                    data[:LL] = data[:LL][2:end]
                    data[:bins] = data[:bins][2:end]
                    if ~isempty(LL)
                        data, cdfs_ii_2, chamfer_ii = resample_data( RR_ii, data )
                    end
                    data[:LL] = LL
                    data[:bins] = bins[ii]

                    cdfs_ii = hcat( cdfs_ii, cdfs_ii_2 )
                else
                    data, cdfs_ii, chamfer_ii = resample_data( RR_ii, data )
                end

                push!( ecdfs, cdfs_ii )
                push!( chamfer_dists, chamfer_ii )
            end
            ecdfs = vcat( ecdfs... )
            chamfer_dists = vcat( chamfer_dists... )
            summary_stats = hcat( ecdfs, chamfer_dists )

            data[:bins] = bins  # Reset bins in data dict
            data[:res_dim][1] = data[:case_dim][2]  # Set resampling dimension correctly
        end
    else
        # TODO
        # Make bins by epochs
        # data, cdf = pairwise_distL(data[:R0], data)
    end

    if data[:log] == "log"
        summary_stats = log.( summary_stats )
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

    data[:muu_data] = muu
    data[:C] = C
    println(size(summary_stats))
    return data, summary_stats
end

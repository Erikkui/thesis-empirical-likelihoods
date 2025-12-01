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
    create_kdtree = kn_e > 0 || data[:chamfer] == 1 ? true : false
    data[:create_kdtree] = create_kdtree
    data[:kn] = kn

    # Initialize bin storage
    data[:bins] = Dict{Symbol, Any}()
    if data[:eCDF] == 1
        data[:bins] = Vector{ Vector{ Vector{Float64} } }(undef, 1 + length(data[:diff_order])*data[:use_diff] )
    end

    if data[:resample] == 1     # Make bins and CDF/chamfer dist by resampling
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
            for ii in eachindex( data[:R0_all] )
                data, _, bins = resample_data( data[:R0_all][ii], data )
                data[:bins][ii] = bins
            end
            data[:bins_done] = true

            LL = data[:LL]  # Store original LL
            bins = data[:bins]  # Store bins for later use
            ecdfs = Vector{ Float64 }(undef, 0)
            chamfer_dists = Vector{ Float64 }(undef, 0)
            for ii in eachindex( data[:R0_all] )
                RR_ii = data[:R0_all][ii]
                ndata_R_ii = size( RR_ii, 2 )
                data[:bins] = bins[ii]

                x_ind = rand( 1:ndata_R_ii, ndata_R_ii )
                unique!( x_ind )
                y_ind = setdiff( 1:ndata_R_ii, x_ind )
                x = @view RR_ii[ :, x_ind ]
                y = @view RR_ii[ :, y_ind ]

                if -1 in data[:LL]          # Signal feature
                    cdfs_ii, _ = empcdf( vec(RR_ii), nx=data[:nbin], x=data[:bins][1] )
                    data[:LL] = data[:LL][2:end]
                    data[:bins] = data[:bins][2:end]
                    if ~isempty(LL)
                        cdfs_ii_2, chamfer_ii = create_summaries( x, y, data, create_kdtree, kn; same = true )
                    end
                    data[:LL] = LL
                    data[:bins] = bins[ii]

                    cdfs_ii = vcat( vec(cdfs_ii), cdfs_ii_2 )

                else
                    cdfs_ii, chamfer_ii = create_summaries( x, y, data, create_kdtree, kn; same = true )
                end

                append!( ecdfs, cdfs_ii )
                append!( chamfer_dists, chamfer_ii )
            end
            data[:muu_data] = vcat( ecdfs, chamfer_dists )
            data[:bins] = bins  # Reset bins in data dict

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

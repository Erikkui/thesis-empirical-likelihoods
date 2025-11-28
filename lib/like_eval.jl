function like_eval( theta, data::Dict{Symbol, Any} )

    # Call the wrapper function to get updated data and simulated CDF
    data[:theta] = theta        # Update data with the current parameter

    if any( theta .<= 0 )
        return ss = 1e10
    end

    data, summary_stats = Wrapper( theta, data )

    # Compute the mean of the summary stats if not done already
    muu_ss = size( summary_stats, 1 ) > 1 ? mean( summary_stats, dims=1 ) : summary_stats

    C = data[:C]
    muu_R0 = data[:muu]

    # Handle numerical issues
    if any(isnan, C)
        error("Covariance matrix contains NaN values.")
    end

    # if cond(C) >= 1e10
        C = C + 1e-6 * I    # Regularization
        if cond(C) >= 1e10
            ss = NaN
            return ss
        end

        delta = vec( muu_ss ) - vec( muu_R0 )

        ss = (delta' * (C \ delta)) + logdet(C)

        return ss
    # else

    # end
end

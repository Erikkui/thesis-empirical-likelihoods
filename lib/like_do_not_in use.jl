function like_do( data )

    np = size( data[:r], 2 )
    cdf = zeros( np, data[:nbin] )

    for i in 1:np
        cdf[i, :], _ = empcdf( data[:r][:, i], nx=data[:nbin], x=data[:bins] )
    end

    if data[:uni] == "log"
        cdf = log.(cdf)
    end

    if nargout == 1
        return cdf
    end

    # # compute muu and C
    # if size(cdf, 1) > 1
    #     muu = vec(mean(cdf, dims=1))
    # else
    #     muu = cdf
    # end

    # if nargout == 2
    #     if data[:C_how] == "don"
    #         C = donsker(muu, data[:N])
    #     elseif data[:C_how] == "cov"
    #         C = cov(cdf)
    #     else
    #         error("Unknown C_how method: $(data[:C_how])")
    #     end

    #     if get(data, :count, 0) == 1
    #         data[:cdf] = muu
    #     end

    #     data[:C] = C
    # end

    # return cdf, data
end

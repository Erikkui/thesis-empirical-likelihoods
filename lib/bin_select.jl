function bin_select( data )
    # Generate bins for empirical cdf calculation

    # Empirical from simulations
    if length(data[:minmax]) == 2
        a = data[:minmax][1]
        b = data[:minmax][2]
    else
        r = data[:r]
        a = minimum(vec(r))
        b = maximum(vec(r))
        nd = length(r)
    end

    delta = (b - a) / 200
    a += delta
    b -= delta

    if data[:uni] == "xax"
        bins = collect( range(a, b, length=data[:nbin]) )
    elseif data[:uni] == "yax"
        nbin_temp = 100
        bins_temp = collect( range(a, b, length=nbin_temp) )

        # Dense ecdf for ecdf
        np = size( data[:r], 2 )
        cdfs = zeros( np, nbin_temp )
        for ii in 1:np
            cdfs[ii, :], _ = empcdf( data[:r][:, ii], nx=nbin_temp, x=bins_temp )
        end

        if size(data[:r], 2) > 1
            mean_cdf = mean( cdfs, dims=1 )
            mean_cdf = vec( mean_cdf )
        else
            mean_cdf = cdfs
        end

        # Inverse CDF for final bins
        bins = invcdf( bins_temp, mean_cdf, data[:nbin], 1)

    elseif data[:uni] == "log"
        R0 = b
        bb = (R0 / a / 1.01)^(1 / data[:nbin])
        bins = R0 .* bb .^ (-data[:nbin]:-1)
    end

    return bins
end

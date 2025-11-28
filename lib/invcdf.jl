function invcdf(x, cdf, nr, cont=1)
    xi = x[:]
    cdfi = cdf[:]
    rr = range(1.01 * minimum(cdf), stop=0.99 * maximum(cdf), length=nr)
    r = zeros( Float64, nr )
    n_xi = length(xi)

    for i in 1:nr

        arg = sum( rr[i] .> cdfi ) + 1
        ind = min( arg, n_xi )

        if ind == 1 || cont == 2
            r[i] = xi[ind]
        else
            r[i] = ( xi[ind] - xi[ind-1] ) / ( cdfi[ind] - cdfi[ind-1] ) *
                   ( rr[i] - cdfi[ind-1] ) + xi[ind-1]
        end
    end

    return r
end

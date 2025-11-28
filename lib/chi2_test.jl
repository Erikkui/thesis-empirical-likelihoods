function chi2_test(Y, muu=nothing, iC=nothing)
    # Chi-square test for normality

    if muu === nothing && iC === nothing  # chi2 normality test for rows of Y
        Yave = mean(Y, dims=1)
        Ystd = std(Y, dims=1)
        Y0 = Y .- Yave
        C = cov(Y)

        # println( "C size: ", size(C))
        # println("Condition number of summary_stats: ", cond(C))
        # println("Minimum singular value of summary_stats: ", minimum(svdvals(C)))
        # println("Maximum singular value of summary_stats: ", maximum(svdvals(C)))
        # println("Rank of summary_stats: ", rank(C))
        # println("Is summary_stats invertible? ", isposdef(C) && rank(summary_stats) == size(C, 1))
        # println("Isposdef: ", isposdef(C))
        # println("rank = size1: ", rank(C), " == ", size(C, 1))


        iC = inv(C)
    else  # test if the rows of Y follow N(muu,C)
        Yave = mean(Y, dims=1)
        Ystd = std(Y, dims=1)
        Y0 = Y .- muu
    end

    nlogl = sum(Y0' .* (iC * Y0'), dims=1)

    # Histogram calculation
    h = fit( Histogram, vec(nlogl), nbins = 50 )
    x = collect( h.edges[1][1:end-1] ) .+ Float64( h.edges[1].step )/2  # Midpoints of histogram bins
    khi = h.weights
    khi_n = khi ./ ( sum(khi) * (x[2] - x[1]) )

    # Chi-square probability function
    chi_pf = pdf( Chisq( size(Y0, 2) ), x )

    return nlogl, iC, x, chi_pf, Yave, Ystd, khi_n
end

function chipf(x, df)
    # CHIPF Chi squared probability density function
    # CHIPF(x,df), x value, df degrees of freedom

    # Note: In Julia, we can directly use the Chi-square distribution
    return pdf.( Chisq(df), x)
end

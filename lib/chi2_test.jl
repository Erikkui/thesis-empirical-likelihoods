function chi2_test(Y, muu=nothing, iC=nothing)
    n, k = size(Y)
    # Chi-square test for normality

    if muu === nothing && iC === nothing  # chi2 normality test for rows of Y
        Yave = mean(Y, dims=1)
        Ystd = std(Y, dims=1)
        Y0 = Y .- Yave
        C = cov(Y)
        iC = inv(C)
    else  # test if the rows of Y follow N(muu,C)
        Yave = mean(Y, dims=1)
        Ystd = std(Y, dims=1)
        Y0 = Y .- muu
    end

    D_sq = sum(Y0' .* (iC * Y0'), dims=1)

    # Theoretical Chi-square PDF
    chisq_dist = Chisq(k)

    # Emprical Histogram (pdf) calculation
    h = fit( Histogram, vec(D_sq), nbins = 50 )
    x = collect( h.edges[1][1:end-1] ) .+ Float64( h.edges[1].step )/2  # Midpoints of histogram bins
    khi = h.weights
    khi_n = khi ./ ( sum(khi) * (x[2] - x[1]) )

    # Chi-square theoretical pdf
    chi_pf = pdf( chisq_dist, x )

    D_sorted = sort( vec(D_sq) )
    probs = (1:n) ./ (n + 1)
    theo_q = quantile.(chisq_dist, probs)

    return x, chi_pf, khi_n, theo_q, D_sorted
end

using Random
using Distributions: Chisq, pdf, quantile
using StatsBase: fit, Histogram
using Statistics: mean, std, cov
using LoopVectorization
using LinearAlgebra: inv
using CairoMakie

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
    h = fit( Histogram, vec(D_sq), nbins = 100 )
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

function empcdf( data::AbstractArray{Float64}; nx=10, x=nothing)

    if isnothing(x)
        m = minimum( data )
        M = maximum( data )
        x = range( m, stop=M, length=nx )
    end

    ndata = length( data )
    cdf_emp = zeros( Float64, (1, nx) )

    @turbo for ii in 1:nx, jj in 1:ndata
            cdf_emp[ii] += ( data[jj] <=  x[ii] )
    end

    cdf_emp /= ndata

    return cdf_emp, x
end

function plot_chi( chi2_x, chi_true, chi_test )
    colors = Makie.to_colormap(:seaborn_dark)
    set_theme!(palette = (; color = colors), fontsize = 24)

    # 3. Initialize Figure
    # Width 1200 to match the side-by-side style of previous plots
    fig = Figure(size = (1100, 400), figure_padding = (5, 5, 5, 10))

    # --- Left Plot: Chi-Square Test ---
    ax_chi = Axis(fig[1, 1],
        # xlabel = L"\mathrm{Squared\ Mahalanobis\ Distance}\ (D^2)",
        # ylabel = L"\mathrm{Density}",
        xgridvisible = true, ygridvisible = true,
        # xticklabelsize = 18, yticklabelsize = 18
    )

    # Plot lines as requested (using explicit colors from your snippet)
    # Slightly thicker lines (2.5) for thesis readability
    lines!(ax_chi, chi2_x, chi_test, color = colors[2], linewidth = 4)
    lines!(ax_chi, chi2_x, chi_true, color = colors[1], linewidth = 4)
        #    label = L"\mathrm{Empirical\ PDF}")

    return fig
end

function main()
    n_cdf = 10000
    ndata = 200
    nbin = 10

    R0 = rand( ndata )
    xmin = 1.05*minimum(R0)
    xmax = 0.95*maximum(R0)
    bins = collect( range( xmin, xmax, length=nbin ) )

    cdfs = zeros( n_cdf, nbin )
    for ii in 1:n_cdf
        x_ind = rand( 1:ndata, ndata )
        unique!( x_ind )

        x_ii = @view R0[ x_ind ]

        cdf_ii, _ = empcdf( x_ii; nx=nbin, x=bins )
        cdfs[ii, :] = vec( cdf_ii )
    end

    x, chi_pf, khi_n, theo_q, D_sorted = chi2_test( cdfs )

    fig = plot_chi(x, chi_pf, khi_n)
    display( fig )
    save("/home/eki/GitHub/thesis-empirical-likelihoods/thesis_figures/example_realisations/example_chi2.pdf", fig)
end

main()

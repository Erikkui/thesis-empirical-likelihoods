function Wrapper( theta, data )

    Ndata = data[:Ndata]
    nsim = data[:nsim]
    nbin = data[:nbin]
    bins = data[:bins]
    R0 = data[:R_all][1, :]
    LL = copy( data[:LL] )

    summary_stats1 = zeros( nsim, nbin )
    summary_stats2 = Matrix{Float64}(undef, nsim, length(LL)*nbin )
    feat_ii = Vector{Float64}()
    for sim in 1:nsim
        R = randn( 1, Ndata )*theta[2] .+ theta[1]
        ecdf, _ = empcdf( R, nx=nbin, x=bins[1] )

        if !isempty(data[:LL])
            D = pairwise(Euclidean(), reshape(R0, 1, :), reshape(R, 1, :); dims=2)
            cdfs = CIL_ID( D, LL )
            for feature in eachindex( cdfs )

                cdf_feature, _ = empcdf( cdfs[feature], nx=nbin, x=bins[1+feature] )
                feat_ii = vcat( feat_ii, vec( cdf_feature ) )

            end

            summary_stats2[sim, :] = feat_ii
        end

        summary_stats1[sim, :] = vec( ecdf )
    end

    # fig = Figure()
    # ax = Axis(fig[1, 1], title="Summary Statistics", xlabel="Bin", ylabel="Value")
    # for sim in 1:nsim
    #     lines!(ax, bins, summary_stats[sim, :])
    # end
    # display(fig)
    if !isempty(data[:LL])
        summary_stats = hcat( summary_stats1, summary_stats2 )
    else
        summary_stats = summary_stats1
    end

    data[:C] = cov( summary_stats )
    return data, summary_stats
end

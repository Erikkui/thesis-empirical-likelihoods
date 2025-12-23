
function OU_main( data, options, model; datafile = nothing )

    if data[:synthetic_data]   # Simulate nepo data sets
        data = create_synth_data( data )
    else
        # To be added?
    end

    data[:data_dim] = size( data[:R0][1], 1 )

    # Generate run identifier
    identifier = generate_identifier( CASE_NAME, data )
    println( "Starting run ", identifier )

    # Create the likelihood
    data, summary_stats = create_likelihood( data )

    # Test gaussianity by khi^2 test
    fig = Figure();
    ax = Axis(fig[1, 1], title="Khi2 test")
    ax2 = Axis(fig[1, 2], title="Q-Q plot")
    # the cdf vectors by the khi2 test:
    nlogl,iC,x,chi_pf,Yave,Ystd,khi_n,theo_q,D_sorted = chi2_test( summary_stats )
    lines!( ax, x, chi_pf, color=:red, label="Chi-square PDF" )
    lines!( ax, x, khi_n, color=:green, label="Chi-square histogram" )
    scatter!(ax2, theo_q, D_sorted, color=:blue, markersize=8)
    ablines!(ax2, 0, 1, color=:red, linestyle=:dash) # 1:1 line
    axislegend(ax; position=:rt)
    display( fig )
    data[ :chi2_x ] = x
    data[ :chi2_pdf ] = chi_pf
    data[ :chi2_hist ] = khi_n
    data[ :qq_theor_q ] = theo_q
    data[ :qq_D_sorted ] = D_sorted
    sleep(1)

    # Evaluate likelihood before MCMC
    theta = data[:params]
    ss = like_eval(theta, data)
    println( "Initial likelihood value: ", ss )

    data[:params] = theta .* randn( size(theta) )  # Perturb initial parameters for MCMC

    # Run MCMC (uncomment the line you want to use; timev for benchmarking, profview for profiling)
    chain, sschain, results = mcmcrun2(model, data, options)
    # @timev chain, sschain, results = mcmcrun2(model, data, options)    # For timing
    # @profview _, _, _ = mcmcrun2(model, data, options)    # For profiling (comment out when not needed)

    # Create figures and save NetCDF if specified
    if data[:save_figures] || data[:show_figures]
        create_figures(chain, theta, FIG_PATH; show=data[:show_figures], savefig=data[:save_figures], case_name=identifier)
    end
    if data[:save_netcdf]
        save_netcdf( chain, sschain, results, data, NETCDF_PATH, identifier=identifier )
    end

    return chain, sschain, results, data
end

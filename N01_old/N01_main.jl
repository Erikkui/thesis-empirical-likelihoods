
function N01_main( data, options, model; datafile = nothing )

    data = create_synth_data( data )

    data[:data_dim] = size( data[:R0][1], 1 )


    # # Test gaussianity by khi^2 test
    # data, summary_stats = create_likelihood( data )
    # fig = Figure();
    # ax = Axis(fig[1, 1], title="Khi2 test")
    # the cdf vectors by the khi2 test:
    # nlogl,iC,x,chi_pf,Yave,Ystd,khi_n = chi2_test( summary_stats )
    # lines!( ax, x, chi_pf, color=:red, label="Chi-square PDF" )
    # lines!( ax, x, khi_n, color=:green, label="Chi-square histogram" )
    # axislegend(ax; position=:rt)
    # display( fig )


    # Generate run identifier
    identifier = generate_identifier( CASE_NAME, data )
    println( "Starting run ", identifier )

    data[:data_dim] = size( data[:R0][1], 1 )  # Number of dimensions in the data

    # Create the likelihood
    data, _ = create_likelihood( data )


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

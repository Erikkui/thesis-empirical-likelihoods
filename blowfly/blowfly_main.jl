
function blowfly_main( data, options, model; datafile = nothing )

    if data[:chamfer] == 0 && data[:eCDF] == 0
        error( "No summary statistics to compute. Set at least either chamfer or eCDF to 1." )
    end

    # Data generation or loading
    if data[:synthetic_data]   # Simulate nepo data sets
        data = create_synth_data( data )
    else
        data = load_data( datafile, data; dataset_ind = data[:dataset_num] )
        data[:synth_init] = data[:R0][1]  # Initial condition from the first time point
        data[:R0] = [ data[:R0] ]

        if data[:use_diff] == 1
            R0_diffs = Vector{ Vector{Matrix{Float64}} }(undef, length(data[:diff_order]) )
            for ii in eachindex(data[:diff_order])
                R0_diffs[ii] = Vector{ Matrix{ Float64 } }(undef, 0)
            end
            R0_diff = calculate_diffs( data[:R0][1], data[:diff_order], data[:synth_dt] )
            push!.( R0_diffs, R0_diff )
            data[:R0_diff] = R0_diffs
        end
    end

    # Generate run identifier
    identifier = generate_identifier( CASE_NAME, data )
    println( "Starting run ", identifier )

    data[:data_dim] = size( data[:R0][1], 1 )  # Number of dimensions in the data

    # Create the bins (BSL) or likelihood and bins (GSL)
    data, _ = create_likelihood( data )

    # Evaluate likelihood before MCMC
    theta = data[:params]
    ss = like_eval(theta, data)
    println( "Initial likelihood value: ", ss )

    # Run MCMC (uncomment the line you want to use; timev for benchmarking, profview for profiling)
    # @profview _, _, _ = mcmcrun2(model, data, options)    # For profiling (comment out when not needed)
    chain, sschain, results = mcmcrun2(model, data, options)
    # @timev chain, sschain, results = mcmcrun2(model, data, options)    # For timing

    # Create figures and save NetCDF if specified
    if data[:save_figures] || data[:show_figures]
        create_figures(chain, theta, FIG_PATH; show=data[:show_figures], savefig=data[:save_figures], case_name=identifier)
    end
    if data[:save_netcdf]
        save_netcdf( chain, sschain, results, data, NETCDF_PATH, identifier=identifier )
    end

    return chain, sschain, results, data
end

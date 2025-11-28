
function save_netcdf( chain,
                      sschain,
                      results,
                      data,
                      save_path;
                      file_ind = 0,
                      identifier = nothing )
    # Save the chain and sschain to a NetCDF file
    # data is a dictionary containing metadata
    # chain: MCMC chain (Nsim x npar matrix)
    # sschain: summary statistics chain (Nsim vector)

    # Path to save the file
    if isnothing( identifier )
        fname = string( CASE_NAME, "_", file_ind, ".nc" )  # Name of the NetCDF file
    else
        fname = string( identifier, ".nc" )  # Name of the NetCDF file
    end

    netcdf_file = joinpath( save_path, fname )

    # Ensure no existing file with the same name, increment file_ind if necessary
    while true
        if isfile( netcdf_file )
            fname = string( identifier, "_", file_ind+1, ".nc" )  # Name of the NetCDF file
            netcdf_file = joinpath( save_path, fname )
            file_ind += 1
        else
            break
        end
    end

    Nsim, npar = size( chain )

    println(netcdf_file)
    ds = Dataset( netcdf_file, "c" )  # "c" = create new file

    # Define dimensions
    defDim( ds, "sample", Nsim )
    defDim( ds, "param", npar )

    # Define variables
    var_chain = defVar( ds, "chain", Float64, ("sample", "param") )
    var_sschain = defVar( ds, "sschain", Float64, ("sample", ) )

    # Assign data
    var_chain[:, :] = chain
    var_sschain[:] = vec( sschain )  # flatten sschain to 1D

    ds.attrib["title"] = "MCMC data for " * CASE_NAME * " parameter estimation"
    # Add metadata from the dictionary "data"
    for (k, v) in data
        try
            if length(v) < 10 && string(k) != "R0"  # Avoid large data like R0
                ds.attrib[ string(k) ] = string(v)  # attributes must be scalar, use string() for safety
            end
        catch e
            continue
        end
    end
    ds.attrib[ "accept" ] = results[:accept]

    close(ds)

    println( "\nNetCDF file saved to: ", netcdf_file )

end

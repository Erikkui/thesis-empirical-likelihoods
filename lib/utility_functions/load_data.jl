function load_data( file_name::String, data::Dict{Symbol, Any}; dataset_ind = 1::Int )
    data_path = joinpath( pwd(), CASE_NAME, "materials", file_name)
    file = CSV.File( data_path )

    dataset_mask = file.set .== dataset_ind
    xdata = file.day[ dataset_mask ]
    ydata = float.( file.count[ dataset_mask ] )
    dx = xdata[2] - xdata[1]

    if size( xdata, 1 ) > size( xdata, 2)
        ydim = size( ydata, 2)
        ydata = reshape( ydata, ydim, : )
    end

    t = length( xdata )

    println( size(ydata))

    data[:R0] = ydata
    data[:R0_full] = ydata
    data[:nobs] = length( ydata )
    data[:synth_N] = t
    data[:synth_init] = ydata[1]
    data[:synth_dt] = dx

    return data
end



# function load_data_test( file_name::String, dataset_ind = 1 )
#     data_path = joinpath( pwd(), "blowfly", "materials", file_name)
#     file = CSV.File( data_path  )

#     dataset_mask = file.set .== dataset_ind
#     xdata = file.day[ dataset_mask ]
#     ydata = file.count[ dataset_mask ]


#     return xdata, ydata
# end


## Testing
# using CSV
# using CairoMakie
# file_name = "blowflies.csv"

# xdata, ydata = load_data_test( file_name, 1 )
# println( length(xdata), "    ", length(ydata) )

# fig = Figure( size = (1280, 600) )
# ax = Axis( fig[1, 1], xlabel = "Days", ylabel = "Count",
#     xlabelsize = 30, ylabelsize = 30, xticklabelsize = 25, yticklabelsize = 25 )
# lines!(ax, xdata, ydata, label = "Blowflies data", color = :blue, linewidth = 2 )

# display(fig)

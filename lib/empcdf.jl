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

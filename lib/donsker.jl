function donsker( mea, N )
    len = length(mea)
    C = Matrix{Float64}(undef, len, len)

    for ii in 1:len
        for jj in 1:len
            C[ii, jj] = min(mea[ii], mea[jj]) - mea[ii] * mea[jj]
        end
    end

    C /= N
    return C
end

function pairwise_turbo(x::AbstractArray, y::AbstractArray)
    n_x = size(x, 2)
    n_y = size(y, 2)
    d = size(x, 1)

    # Pre-allocate output
    out = Vector{Float64}(undef, n_x*n_y)

    # @tturbo applies threading and SIMD
    @tturbo for j in 1:n_y
        for i in 1:n_x
            s = 0.0
            for k in 1:d
                diff = x[k, i] - y[k, j]
                s += diff * diff
            end
            out[ (j-1)*n_y + i ] = s
        end
    end
    return out
end

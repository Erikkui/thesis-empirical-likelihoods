function calculate_diffs( R, diff_order, dt)

    ndim, ndata = size(R)
    R0_diffs = Vector{ Matrix{Float64} }(undef, length(diff_order) )
    ord_max = maximum( diff_order )
    R0_diff_ind = 1

    for dim in 1:ord_max
        # Central difference: (f(x+1) - f(x-1))/2
        data_prev = @view R[ :, 1:end-2 ]
        data_next = @view R[ :, 3:end ]
        diff_temp = (data_next - data_prev) ./ (2*dt)

        R = diff_temp

        if dim in diff_order
            R0_diffs[ R0_diff_ind ] = diff_temp
            R0_diff_ind += 1
        end

    end

    return R0_diffs
end

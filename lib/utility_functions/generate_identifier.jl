function generate_identifier( model::String, data )

    case = data[:case]
    C_how = data[:C_how]
    nobs = data[:nobs]
    nepo = data[:nepo]
    nsim, nres = data[:case_dim]

    identifier = string( model,
                         "_", case,
                         "_", C_how,
                         "_n", nobs,
                         "n", nepo,
                         "_s", nsim,
                         "_r", nres )

    if data[:use_diff] == 1
        diff_order = join( string.( data[:diff_order] ) )
        identifier *= "_dt"*diff_order
    end

    if data[:eCDF] == 1
        LL = join( string.(data[:LL]) )
        ecdf = "_eCDF"*LL
        identifier *= ecdf
    end

    if data[:chamfer] == 1
        chamfer_k = join( string.(data[:chamfer_k]) )
        chamf = "_chamfer"*chamfer_k
        identifier *= chamf
    end

    if data[:synthetic_data]
        identifier *= "_synth"
    else
        setnum  = string( data[:dataset_num] )
        identifier *= "_set"*setnum
    end

    return identifier
end

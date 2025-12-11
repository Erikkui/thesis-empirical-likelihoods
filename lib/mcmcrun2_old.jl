
function mcmcrun2_old(model::Dict{Symbol, Any},
                  data::Dict{Symbol, Any},
                  options::Dict{Symbol, Any},
                  results_prev=nothing)

    chain_length = options[:nsimu]
    adaptint = options[:adapt_int]
    update_int = options[:update_int]
    sigma2 = model[:sigma2]
    ssfun = model[:ssfun]
    parameters = data[:theta]

    if isnothing( results_prev )
        qcov = options[:qcov]
    else        # continue from a previous run, using the last values
        parameters = copy( results_prev[:last] )
        qcov = copy( results_prev[:qcov] )
        Chain = copy( results_prev[:Chain] )
    end

    R = cholesky( qcov ).L
    n_params = length( parameters )
    params_old = copy( parameters )
    n_rejected = 0  # initialize count for rejections
    n_rejected_stuck = 0

    chain = zeros( chain_length, n_params )
    sschain = zeros( chain_length )

    chain[1, :] = params_old

    ss = ssfun( params_old, data )  # first SS value


    sschain[1] = ss
    status = 1

    iter = ProgressBar( 2:chain_length )
    ii = 1
    # try
        for ii in iter  # Simulation loop

            rr = randn( Float64, n_params )
            param_new = params_old + R * rr # New parameter candidate

            ss_old = ss #+ logprior_old  # old SS
            ss_new = ssfun( param_new, data ) #+ logprior_proposal  # new SS

            random_accept = exp( -0.5*(ss_new - ss_old) / sigma2 )

            if ss_new < ss_old || rand() < random_accept  # Accept proposal?
                chain[ ii, : ] = param_new  # accept
                params_old = copy( param_new )
                ss = ss_new
                status = 1  # for accept
            else
                if status == 1
                    n_rejected_stuck = 0
                end  # start counting consequent rejections
                status = 0
                n_rejected_stuck += 1
                if n_rejected_stuck > update_int
                    ss_new = ssfun( params_old, data )  # update value at params_old

                    ss = ss_new
                    chain[ ii, : ] = params_old
                    n_rejected_stuck = 0
                else
                    chain[ ii, : ] = params_old  # reject
                    n_rejected += 1
                end
            end

            # ADAPT proposal covariance
            if ii > 200  # no adaptation too early
                if ii % adaptint == 0
                    if !isnothing( results_prev )
                        C = vcat( Chain, chain[1:ii, :] )
                    else
                        C = chain[ 1:ii, : ]
                    end
                    qcov = cov(C)
                    if cond( qcov ) < 1e+10
                        R = cholesky( qcov ).L
                    end
                end
            end

            sschain[ii] = ss

            set_description( iter, "Running MCMC:" )
        end
    # catch e
    #     println( "Error occurred during MCMC:: ", e, "\nReturning results up to the last successful iteration." )
    #     accept = 1 - n_rejected / chain_length  # acceptance rate
    #     results_new = Dict(
    #     :accept => accept,
    #     :last => vec( chain[end, :] ),
    #     :qcov => qcov
    #     )

    #     if !isnothing( results_prev )
    #         results_new[:Chain] = vcat( Chain, chain )
    #     else
    #         results_new[:Chain] = chain
    #     end

    #     return chain, sschain, results_new
    # end

    accept = 1 - n_rejected / chain_length  # acceptance rate
    results_new = Dict(
        :accept => accept,
        :last => vec( chain[end, :] ),
        :qcov => qcov
    )

    if !isnothing( results_prev )
        results_new[:Chain] = vcat( Chain, chain )
    else
        results_new[:Chain] = chain
    end

    return chain, sschain, results_new
end

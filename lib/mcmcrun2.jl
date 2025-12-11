function mcmcrun_stochastic(model::Dict{Symbol, Any},
                            data::Dict{Symbol, Any},
                            options::Dict{Symbol, Any},
                            results_prev=nothing)

    # --- Setup ---
    chain_length = options[:nsimu]
    adaptint = options[:adapt_int]
    update_int = options[:update_int] # Threshold for re-evaluation
    sigma2 = model[:sigma2]
    ssfun = model[:ssfun]

    # --- Initialization ---
    if isnothing(results_prev)
        parameters = data[:theta]
        qcov = options[:qcov]

        # Online statistics init
        global_mean = copy(parameters)
        global_cov = copy(qcov)
        n_history = 0
    else
        parameters = copy(results_prev[:last])
        qcov = copy(results_prev[:qcov])

        # Try to recover history, or reset if missing (safer for adaptation)
        global_mean = get(results_prev, :mean, copy(parameters))
        global_cov = get(results_prev, :cov, copy(qcov))
        n_history = get(results_prev, :n_history, 0)
    end

    n_params = length(parameters)

    # Pre-allocate
    chain = zeros(Float64, chain_length, n_params)
    sschain = zeros(Float64, chain_length)

    params_current = copy(parameters)
    ss_current = ssfun( params_current, data )

    # Proposal distribution setup
    # Optimal scaling: (2.4^2 / d)
    sd = (2.4^2) / n_params
    epsilon = 1e-6

    # Initial Cholesky
    R = cholesky(qcov).L

    n_accepted = 0
    n_stuck    = 0

    iter = ProgressBar(1:chain_length)

    for ii in iter

        # 1. Propose
        noise = randn(n_params)
        params_proposal = params_current + R * noise

        ss_proposal = ssfun(params_proposal, data)

        # 2. Metropolis accept/reject
        log_ratio = -0.5 * (ss_proposal - ss_current) / sigma2

        accepted = false
        if log_ratio >= 0 || log( rand() ) < log_ratio
            accepted = true
        end

        if accepted
            params_current = params_proposal
            ss_current = ss_proposal
            n_accepted += 1
            n_stuck = 0 # Reset stuck counter on move
        else
            # Rejected
            n_stuck += 1
            if n_stuck > update_int
                # If mcmc gets stuck in "too good" proposal
                ss_recalc = ssfun(params_current, data)
                ss_current = ss_recalc
                n_stuck = 0 # Reset counter
            end
        end

        # Store
        chain[ ii, : ] = params_current
        sschain[ii] = ss_current

        # 3. Adaptive Update (Recursive Welford)
        n_history += 1
        delta = params_current - global_mean
        global_mean .+= delta ./ n_history
        delta2 = params_current - global_mean

        if n_history > 1
             factor1 = (n_history - 2) / (n_history - 1)
             factor2 = 1.0 / (n_history - 1)
             @. global_cov = factor1 * global_cov + factor2 * (delta * delta2')
        end

        # 4. Update Proposal Cholesky
        if ii > 100 && ii % adaptint == 0
            prop_cov = sd .* global_cov
            for k in 1:n_params
                prop_cov[k,k] += epsilon
            end
            try
                R = cholesky( Symmetric(prop_cov) ).L
            catch
                # Ignore non-posdef failures temporarily
            end
        end

        set_description(iter, "Acc: $(round(n_accepted/i, digits=2)) SS: $(round(ss_current, digits=2))")
    end

    # Save results
    results_new = Dict(
        :accept    => n_accepted / chain_length,
        :last      => chain[end, :],
        :qcov      => global_cov,
        :mean      => global_mean,
        :n_history => n_history,
        :Chain     => chain
    )

    if !isnothing(results_prev)
        results_new[:Chain] = vcat(results_prev[:Chain], chain)
    end

    return chain, sschain, results_new
end

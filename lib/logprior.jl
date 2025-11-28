function logprior( theta )

    # Define the prior bounds for each parameter
    priors = [ (:uniform, 0.02, 1.0),
              (:uniform, 3, 30),
              (:uniform, 10, 1000),
              (:uniform, 0.01, 0.5),
              (:normal, 14, 5),
              (:uniform, 0.01, 0.5) ]

    logprior_value = 0.0
    for (ii, prior) in enumerate(priors)
        prior_type, a, b = prior
        if prior_type == :uniform
            if theta[ii] < a || theta[ii] > b
                # println("unif reject")
                return -Inf  # outside prior support
            end
        elseif prior_type == :normal
            logprior_value += logpdf( Normal(a, b), theta[ii] )
            # println("normal prior: ", logprior_value)
        end
    end

    return logprior_value
end

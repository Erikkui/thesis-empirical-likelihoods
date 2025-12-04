function lorenz3!(du, u, params, t)
    du[1] = params[1] * (u[2] - u[1])
    du[2] = u[1] * (params[2] - u[3]) - u[2]
    du[3] = u[1] * u[2] - params[3] * u[3]
    nothing
end

function lorenz_solve( init::Vector{Float64}, theta::Vector{Float64}, Ndata; dt = 1 )
    t = (0:Ndata)*dt
    s0 = init .* (1 .+ 0.01 .* randn(3))
    # Define ODE problem and solve
    prob = ODEProblem( lorenz3!, s0, (t[1], t[end]), theta)
    sol = solve(prob, Tsit5(), saveat=t)

    # Drop first 10 points
    R0 = Array( sol )[ :, 11:end ]

    return R0
end

function lorenz3!(du, u, params, t)
    du[1] = params[1] * (u[2] - u[1])
    du[2] = u[1] * (params[2] - u[3]) - u[2]
    du[3] = u[1] * u[2] - params[3] * u[3]
    nothing
end

function lorenz_static(u, p, t)
    σ, ρ, β = p
    x, y, z = u

    dx = σ * (y - x)
    dy = x * (ρ - z) - y
    dz = x * y - β * z

    return SVector(dx, dy, dz)
end

function lorenz_solve(init, theta, Ndata; dt=1.0)
    # 2. Use StaticArrays for Inputs
    # Generate noise element-wise to avoid allocating a Vector{Float64}
    noise = SVector(randn(), randn(), randn())
    u0 = SVector{3}(init) .* (1.0 .+ 0.01 .* noise)
    p  = SVector{3}(theta)

    # 3. Smart Saving (No Slicing)
    # Instead of solving everything and slicing later,
    # use `saveat` to only store the points you strictly need.
    # User drops first 10 steps (indices 1:10). Start saving at t = 10*dt.
    t_start = 10.0 * dt
    t_end   = Ndata * dt

    prob = ODEProblem(lorenz_static, u0, (0.0, t_end), p)

    # save_start=false ensures we don't accidentally save t=0
    sol = solve(prob, Tsit5(); saveat=t_start:dt:t_end, save_start=false)

    # 4. Pipeline efficiency
    # sol.u is ALREADY a Vector{SVector{3, Float64}}.
    # This matches the input format your distance function needs.
    # If you absolutely need a Matrix, use stack(sol.u).
    return stack(sol.u)
end

function lorenz_solve_2( init::Vector{Float64}, theta::Vector{Float64}, Ndata; dt = 1 )
    t = (0:Ndata)*dt
    s0 = init .* (1 .+ 0.01 .* randn(3))
    # Define ODE problem and solve
    prob = ODEProblem( lorenz3!, s0, (t[1], t[end]), theta)
    sol = solve(prob, Tsit5(), saveat=t)

    # Drop first 10 points
    R0 = Array( sol )[ :, 11:end ]

    return R0
end

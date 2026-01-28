using OrdinaryDiffEq

function simulate_observer(sys::LinearSystem, K::Vector{T};
                          x0_true::Vector{T},
                          x0_obs::Vector{T},
                          tspan::Tuple{Float64, Float64} = (0.0, 10.0),
                          u::Function = t -> zeros(T, sys.n),
                          B::Union{Matrix{T}, Nothing} = nothing,
                          alg = Tsit5(),
                          kwargs...) where T<:Real
    
    A = sys.A
    C = sys.C
    n = sys.n
    
    if length(x0_true) != n
        error("x0_true must have length $n, got $(length(x0_true))")
    end
    if length(x0_obs) != n
        error("x0_obs must have length $n, got $(length(x0_obs))")
    end
    if length(K) != n
        error("K must have length $n, got $(length(K))")
    end
    
    # Combined state: [x_true; x_obs]
    z0 = vcat(x0_true, x0_obs)
    
    # Combined dynamics
    function observer_dynamics!(dz, z, p, t)
        # Extract states
        x_true = z[1:n]
        x_obs = z[n+1:2*n]
        
        # True system measurement
        y = dot(C, x_true)
        
        # Innovation (measurement residual)
        y_obs = dot(C, x_obs)
        innovation = y - y_obs
        
        # Input
        u_t = u(t)
        
        # True system dynamics: dx/dt = A*x + B*u
        if B !== nothing
            dz[1:n] = A * x_true + B * u_t
        else
            dz[1:n] = A * x_true
        end
        
        # Observer dynamics: dx̂/dt = A*x̂ + B*u + K*(y - C*x̂)
        if B !== nothing
            dz[n+1:2*n] = A * x_obs + B * u_t + K * innovation
        else
            dz[n+1:2*n] = A * x_obs + K * innovation
        end
    end
    
    # Create and solve ODE problem
    ode_prob = ODEProblem(observer_dynamics!, z0, tspan)
    sol = solve(ode_prob, alg; kwargs...)
    
    # Extract results
    t = sol.t
    x_true = reduce(hcat, [sol.u[i][1:n] for i in 1:length(sol)])
    x_obs = reduce(hcat, [sol.u[i][n+1:2*n] for i in 1:length(sol)])
    e = x_true - x_obs
    
    # Compute measurements
    y = [dot(C, x_true[:, i]) for i in 1:size(x_true, 2)]
    
    # Observer dynamics matrix
    C_row = reshape(C, 1, n)
    A_obs = A + K * C_row
    
    return (
        sol = sol,
        t = t,
        x_true = x_true,
        x_obs = x_obs,
        e = e,
        y = y,
        A_obs = A_obs,
        K = K
    )
end


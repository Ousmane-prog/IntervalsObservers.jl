function evaluate_f(f_vec, t, x, y, n)
    f_eval = similar(x)
    @inbounds for i in 1:n 
        f_eval[i] = f_vec[i](t, x, y)
    end 
    return f_eval  
end


function system_dynamics!(dx, x, p, t)
    sys, y_func = p

    y = y_func(t)
    f_val = evaluate_f(sys.f, t, x, y, sys.n)
    dx .= sys.A * x + f_val
end

function non_linear_interval_observer!(dz, z, p, t)
    sys, K, y_func = p 
    n = sys.n

    x_plus = view(z, 1:n)
    x_minus = view(z, n+1:2n)

    y = y_func(t)
    f_plus_val = evaluate_f(sys.f_plus, t, x_plus, y, n)
    f_minus_val = evaluate_f(sys.f_minus, t, x_minus, y, n)
    dz[1:n] .= sys.A * x_plus .+ f_plus_val + K * (y - dot(sys.C, x_plus))
    dz[n+1:2n] .= sys.A * x_minus .+ f_minus_val + K * (y - dot(sys.C, x_minus))
end
"""
    get_state(sol, n::Int)

Extract the true state trajectory `x(t)` from the interval observer solution.

# Arguments
- `sol`: ODE solution from interval observer simulation
- `n::Int`: System dimension

# Returns
- `Matrix`: State trajectory (n×T) where T is the number of time points
"""
function get_state(sol, n::Int)
    return sol[1:n, :]
end

"""
    get_lower(sol, n::Int)

Extract the lower bound trajectory `x⁻(t)` from the interval observer solution.

# Arguments
- `sol`: ODE solution from interval observer simulation
- `n::Int`: System dimension

# Returns
- `Matrix`: Lower bound trajectory (n×T) where T is the number of time points
"""
function get_upper(sol, n::Int)
    return sol[n+1:2n, :]
end

"""
    get_upper(sol, n::Int)

Extract the upper bound trajectory `x⁺(t)` from the interval observer solution.

# Arguments
- `sol`: ODE solution from interval observer simulation
- `n::Int`: System dimension

# Returns
- `Matrix`: Upper bound trajectory (n×T) where T is the number of time points
"""
function get_lower(sol, n::Int)
    return sol[2n+1:3n, :]
end

"""
    get_upper_nonlinear(sol, n::Int)

Extract the upper bound trajectory `x⁺(t)` from a nonlinear interval observer solution.

# Arguments
- `sol`: ODE solution from nonlinear interval observer simulation
- `n::Int`: System dimension

# Returns
- `Matrix`: Upper bound trajectory (n×T) where T is the number of time points
"""
function get_upper_nonlinear(sol, n::Int)
    return sol[1:n, :]
end

"""
    get_lower_nonlinear(sol, n::Int)

Extract the lower bound trajectory `x⁻(t)` from a nonlinear interval observer solution.

# Arguments
- `sol`: ODE solution from nonlinear interval observer simulation
- `n::Int`: System dimension

# Returns
- `Matrix`: Lower bound trajectory (n×T) where T is the number of time points
"""
function get_lower_nonlinear(sol, n::Int)
    return sol[n+1:2n, :]
end


function positive_interval_gain(sys::Union{LinearSystem, NonLinearSystem})
    # sys.positive || error("interval observers require a positive system")

    n = sys.n
    C = sys.C

    K = zeros(n)
    for i in 1:n
        if C[i] > 0
            K[i] = 0.1
            # K[i] = 0
        end 
    end
    return K
end


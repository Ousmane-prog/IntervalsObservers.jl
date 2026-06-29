using DifferentialEquations

"""
    build_nonlinear_interval_problem(obs, x0_plus, x0_minus, tspan; x0=nothing, f_true=nothing)

Construct an ODE problem for a nonlinear interval observer system.

Builds the coupled ODE problem for the observer dynamics:
  ẋxᴶ = A*xᴶ + fᴶ(t,y) + K*(y - C*xᴶ)
  ẋxᴵ = A*xᴵ + fᴵ(t,y) + K*(y - C*xᴵ)

Optionally tracks the true state when x0 and f_true are provided.

# Arguments
- `obs::IntervalObserver`: Observer with system and gain
- `x0_plus::Vector`: Initial upper bound
- `x0_minus::Vector`: Initial lower bound
- `tspan::Tuple`: Time span (t0, tf)
- `x0::Union{Nothing, Vector}`: Initial true state (optional)
- `f_true::Union{Nothing, Vector}`: True nonlinearity functions for each state (optional)

# Returns
- `ODEProblem`: ODE problem ready for solving with DifferentialEquations.jl

# Note
If x0 and f_true are provided, the ODE state includes the true state:
  z = [x; x_plus; x_minus] (length 3n)
Otherwise:
  z = [x_plus; x_minus] (length 2n)
"""
function build_nonlinear_interval_problem(
    obs::IntervalObserver,
    x0_plus::Vector,
    x0_minus::Vector,
    tspan;
    x0::Union{Nothing, Vector} = nothing,
    f_true::Union{Nothing, Vector} = nothing
)
    sys = obs.sys
    n = sys.n
    track_true_state = !isnothing(x0)

    if track_true_state && isnothing(f_true)
        error("Tracking true state requires f_true.")
    end

    if track_true_state
        X0 = vcat(x0, x0_plus, x0_minus)
    else
        X0 = vcat(x0_plus, x0_minus)
    end

    function interval_dynamics!(dX, X, p, t)
        A = sys.A
        C = sys.C
        K = obs.K

        if track_true_state
            x       = @view X[1:n]
            x_plus  = @view X[n+1:2n]
            x_minus = @view X[2n+1:3n]

            y = dot(C, x)
        else
            x_plus  = @view X[1:n]
            x_minus = @view X[n+1:2n]

            y = dot(C, x_plus)
        end

        f_plus_val  = similar(x_plus)
        f_minus_val = similar(x_minus)

        @inbounds for i in 1:n
            f_plus_val[i]  = sys.f_plus[i](t, y)
            f_minus_val[i] = sys.f_minus[i](t, y)
        end

        if track_true_state
            f_true_val = similar(x)

            @inbounds for i in 1:n
                f_true_val[i] = f_true[i](t, y)
            end

            # actual system
            dX[1:n] .= A*x + f_true_val

            # upper observer
            dX[n+1:2n] .= A*x_plus +
                          f_plus_val +
                          K*(y - dot(C, x_plus))

            # lower observer
            dX[2n+1:3n] .= A*x_minus +
                           f_minus_val +
                           K*(y - dot(C, x_minus))
        else
            dX[1:n] .= A*x_plus +
                       f_plus_val +
                       K*(y - dot(C, x_plus))

            dX[n+1:2n] .= A*x_minus +
                          f_minus_val +
                          K*(y - dot(C, x_minus))
        end
    end

    return ODEProblem(interval_dynamics!, X0, tspan)
end
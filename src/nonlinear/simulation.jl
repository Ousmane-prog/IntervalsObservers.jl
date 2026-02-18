function build_nonlinear_interval_problem(
    obs::IntervalObserver,
    x0_plus::Vector,
    x0_minus::Vector,
    tspan;
    x0::Union{Nothing, Vector} = nothing
)

    sys = obs.sys
    n = sys.n
    track_true_state = !isnothing(x0)

    # Build initial condition based on whether x0 is provided
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
            x = @view X[1:n]
            x_plus = @view X[n+1:2n]
            x_minus = @view X[2n+1:3n]
            
            # Measured output from the true system
            y = dot(C, x)
        else
            x_plus = @view X[1:n]
            x_minus = @view X[n+1:2n]
            
            # Measured output from observer bounds (use upper bound as measurement)
            y = dot(C, x_plus)
        end

        # Evaluate nonlinear vector functions
        f_plus_val = [obs.f_plus[i](t, y) for i in 1:n]
        f_minus_val = [obs.f_minus[i](t, y) for i in 1:n]

        if track_true_state

            dX[1:n] .= A*x + f_plus_val 

          
            dX[n+1:2n] .= A*x_plus +
                          f_plus_val +
                          K*(y - dot(C, x_plus))

            
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



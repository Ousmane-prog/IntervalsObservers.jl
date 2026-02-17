function build_nonlinear_interval_problem(
    obs::IntervalObserver,
    x0_plus::Vector,
    x0_minus::Vector,
    tspan
)

    sys = obs.sys
    n = sys.n

    X0 = vcat(x0_plus, x0_minus)

    function interval_dynamics!(dX, X, p, t)

        A = sys.A
        C = sys.C
        K = obs.K

        x_plus = @view X[1:n]
        x_minus = @view X[n+1:2n]

        # Measured output from observer bounds (use upper bound as measurement)
        y = dot(C, x_plus)

        # Evaluate nonlinear vector functions
        f_plus_val = [obs.f_plus[i](t, y) for i in 1:n]
        f_minus_val = [obs.f_minus[i](t, y) for i in 1:n]

        # UPPER OBSERVER
        dX[1:n] .= A*x_plus +
                   f_plus_val +
                   K*(y - dot(C, x_plus))

        # LOWER OBSERVER
        dX[n+1:2n] .= A*x_minus +
                      f_minus_val +
                      K*(y - dot(C, x_minus))
    end

    return ODEProblem(interval_dynamics!, X0, tspan)
end 



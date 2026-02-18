using DifferentialEquations

function solve(
    sys::LinearSystem,
    x0 :: Vector,
    xl0:: Vector,
    xu0::Vector,
    tspan:: Tuple{Real, Real},
    solver = Tsit5()
)
    n = sys.n
    A = sys.A
    C = sys.C
    check_Metzler_Matrix(A)
    validate_initial_bounds(x0, xl0)

    K = positive_interval_gain(sys)
    z0 = vcat(x0, xu0, xl0)
    p = (sys.A, sys.C, K, n)

    prob = ODEProblem(Linear_syst_int_obs_ode!, z0, tspan, p)
    sols = DifferentialEquations.solve(prob, solver)

    return sols
end

function solve(
    obs::IntervalObserver,
    x0_plus::Vector,
    x0_minus::Vector,
    tspan::Tuple{Real, Real};
    x0::Union{Nothing, Vector} = nothing,
    solver = Tsit5()
)
    """
        solve(obs::IntervalObserver, x0_plus, x0_minus, tspan; x0=nothing, solver)
    
    Solve a nonlinear interval observer problem.
    
    # Arguments
    - `obs::IntervalObserver`: Interval observer with system, gain, and bounding functions
    - `x0_plus::Vector`: Initial upper bound on state
    - `x0_minus::Vector`: Initial lower bound on state
    - `tspan::Tuple`: Time span (t0, tf)
    - `x0::Vector` (optional): Initial true state (if available)
    - `solver`: ODE solver (default: Tsit5)
    
    # Returns
    - Solution object containing state trajectories
    """
    validate_initial_bounds(x0_minus, x0_plus)
    
    prob = build_nonlinear_interval_problem(obs, x0_plus, x0_minus, tspan; x0=x0)
    sol = DifferentialEquations.solve(prob, solver)
    
    return sol
end
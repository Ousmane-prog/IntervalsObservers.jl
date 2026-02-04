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
    validate_initial_bounds(x0, xl0, xu0)

    K = positive_interval_gain(sys)
    z0 = vcat(x0, xu0, xl0)
    p = (sys.A, sys.C, K, n)

    prob = ODEProblem(Linear_syst_int_obs_ode!, z0, tspan, p)
    sols = DifferentialEquations.solve(prob, solver)

    return sols
end

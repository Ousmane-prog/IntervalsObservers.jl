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
    sys.positive || error("interval observers require a positive system")
    all(xl0 .<= x0 .<= xu0) || error("Initial conditions must satisfy x⁻₀ ≤ x₀ ≤ x⁺₀")

    K = positive_interval_gain(sys)
    z0 = vcat(x0, xl0, xu0)
    p = (sys.A, sys.C, K, n)

    prob = ODEProblem(Linear_syst_int_obs_ode!, z0, tspan, p)
    sols = DifferentialEquations.solve(prob, solver)

    return sols
end

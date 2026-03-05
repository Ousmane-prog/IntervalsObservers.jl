using DifferentialEquations

function solve(
    sys::LinearSystem,
    x0 :: Vector,
    xl0:: Vector,
    xu0::Vector,
    K::Union{Nothing, Vector},
    tspan:: Tuple{Real, Real},
    solver = Tsit5()
)
    n = sys.n
    A = sys.A
    C = sys.C
    check_Metzler_Matrix(A)
    validate_initial_bounds(x0, xl0)

    z0 = vcat(x0, xu0, xl0)
    p = (sys.A, sys.C, K, n)

    prob = ODEProblem(Linear_syst_int_obs_ode!, z0, tspan, p)
    sols = DifferentialEquations.solve(prob, solver)

    return sols
end

function solve(
    sys::NonLinearSystem,
    K::Vector,
    f_plus::Vector,
    f_minus::Vector,
    x0_plus::Vector,
    x0_minus::Vector,
    tspan::Tuple{Real, Real};
    x0::Union{Nothing, Vector} = nothing,
    solver = Tsit5()
)

    A = sys.A
    C = sys.C
    A_minus_KC = A - K * reshape(C, 1, :)
    validate_initial_bounds(x0_minus, x0_plus)
    if monotone_dynamic(A_minus_KC)

        obs = IntervalObserver(sys, K, f_plus, f_minus)
        prob = build_nonlinear_interval_problem(
            obs, x0_plus, x0_minus, tspan; x0=x0
        )
        prob = build_nonlinear_interval_problem(obs, x0_minus, x0_plus, tspan; x0=x0) 
    else 

        F = eigen(A_minus_KC)
        T = F.vectors
        T_inv = inv(T)
        # x0_plus_new = transform_initial_condition(x0_plus, x0_minus, T_inv)
        z0_minus, z0_plus = transform_interval(T_inv, x0_minus, x0_plus)
        
        x0_new = nothing
        if x0 !== nothing
            z0 = T_inv * x0
        end

        obs = IntervalObserver(sys, K, f_plus, f_minus)
        prob = build_nonlinear_interval_problem(obs, z0_plus, z0_minus, tspan; x0=x0_new)
    end  
    
    sol = DifferentialEquations.solve(prob, solver)

    return sol
end
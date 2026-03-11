# using DifferentialEquations

# function solve(
#     sys::LinearSystem,
#     x0 :: Vector,
#     xl0:: Vector,
#     xu0::Vector,
#     K::Union{Nothing, Vector},
#     tspan:: Tuple{Real, Real},
#     solver = Tsit5()
# )
#     n = sys.n
#     A = sys.A
#     C = sys.C
#     check_Metzler_Matrix(A)
#     validate_initial_bounds(x0, xl0)

#     z0 = vcat(x0, xu0, xl0)
#     p = (sys.A, sys.C, K, n)

#     prob = ODEProblem(Linear_syst_int_obs_ode!, z0, tspan, p)
#     sols = DifferentialEquations.solve(prob, solver)

#     return sols
# end

# function solve(
#     sys::NonLinearSystem,
#     K::Vector,
#     f_plus::Vector,
#     f_minus::Vector,
#     x0_plus::Vector,
#     x0_minus::Vector,
#     tspan::Tuple{Real, Real};
#     x0::Union{Nothing, Vector} = nothing,
#     solver = Tsit5()
# )

#     A = sys.A
#     C = sys.C
#     A_minus_KC = A - K * reshape(C, 1, :)
#     validate_initial_bounds(x0_minus, x0_plus)
#     if monotone_dynamic(A_minus_KC)

#         obs = IntervalObserver(sys, K, f_plus, f_minus)
#         prob = build_nonlinear_interval_problem(
#             obs, x0_plus, x0_minus, tspan; x0=x0
#         )
#         prob = build_nonlinear_interval_problem(obs, x0_minus, x0_plus, tspan; x0=x0) 
#     else 

#         F = eigen(A_minus_KC)
#         T = F.vectors
#         T_inv = inv(T)
#         # x0_plus_new = transform_initial_condition(x0_plus, x0_minus, T_inv)
#         z0_minus, z0_plus = transform_interval(T_inv, x0_minus, x0_plus)
        
#         x0_new = nothing
#         if x0 !== nothing
#             z0 = T_inv * x0
#         end

#         obs = IntervalObserver(sys, K, f_plus, f_minus)
#         prob = build_nonlinear_interval_problem(obs, z0_plus, z0_minus, tspan; x0=x0_new)
#     end  
    
#     sol = DifferentialEquations.solve(prob, solver)

#     return sol
# end

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

# function solve(
#     sys::NonLinearSystem,
#     K::Vector,
#     f_plus::Vector,
#     f_minus::Vector,
#     x0_plus::Vector,
#     x0_minus::Vector,
#     tspan::Tuple{Real, Real};
#     x0::Union{Nothing, Vector} = nothing,
#     solver = Tsit5()
# )

#     A = sys.A
#     C = sys.C
#     A_minus_KC = A - K * reshape(C, 1, :)
#     validate_initial_bounds(x0_minus, x0_plus)
#     if monotone_dynamic(A_minus_KC)

#         obs = IntervalObserver(sys, K, f_plus, f_minus)
#         prob = build_nonlinear_interval_problem(
#             obs, x0_plus, x0_minus, tspan; x0=x0
#         )
#     else 

#         M, M_inv, D = diagonalizing_change_of_basis(A, C, K)
#         # x0_plus_new = transform_initial_condition(x0_plus, x0_minus, T_inv)
#         z0_minus, z0_plus, _z0 = transform_initial_condition(x0_plus, x0_minus, M_inv)
#         println("Transformed initial bounds: z0_minus = ", z0_minus, ", z0_plus = ", z0_plus)
#         x0_new = nothing
#         if x0 !== nothing
#             z0 = M_inv * x0
#         end
#         A_z = M * A * M_inv
#         # C_z = C * T
#         C_z = C * M_inv 
#         f_plus_z = z -> M * f_plus(T * z)
#         f_minus_z = z -> M * f_minus(T * z)

#         new_sys   = NonLinearSystem(A_z, C_z, f_plus_z, f_minus_z; check_metzler=false)
#         obs = IntervalObserver(new_sys, K, f_plus_z, f_minus_z)
#         prob = build_nonlinear_interval_problem(obs, z0_plus, z0_minus, tspan; x0=x0_new)
#     end  
    
#     sol = DifferentialEquations.solve(prob, solver)

#     return sol
# end


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

    validate_initial_bounds(x0_minus, x0_plus)

    A_minus_KC = A - K * reshape(C, 1, :)

    if monotone_dynamic(A_minus_KC)
        obs = IntervalObserver(sys, K, f_plus, f_minus)
        prob = build_nonlinear_interval_problem(
            obs, x0_plus, x0_minus, tspan; x0=x0
        )
    else
        # M diagonalizes A-KC in the sense M(A-KC)M^{-1} = D
        M, M_inv, D = diagonalizing_change_of_basis(A, C, K)

        # Transform initial interval and exact initial condition using z = Mx
        z0_minus, z0_plus = transform_interval(M, x0_minus, x0_plus)
        println("\n⚠️  CHANGE OF BASIS APPLIED:")
        println("   The system was transformed using: z = M*x")
        println("   where M is the change-of-basis matrix:")
        println("   M = ")
        display(M)

        z0 = x0 === nothing ? nothing : M * x0

        # Transform system
        A_z = M * A * M_inv
        C_z = vec((C') * M_inv)

        # Transform nonlinear bounds: f_z(t,y) = M f(t,y)
        f_plus_z  = transform_function_vector(f_plus,  M)
        f_minus_z = transform_function_vector(f_minus, M)

        # Transform observer gain
        K_z = M * K

        new_sys = NonLinearSystem(
            A_z, C_z, f_plus_z, f_minus_z;
            check_metzler = false
        )

        obs = IntervalObserver(new_sys, K_z, f_plus_z, f_minus_z)

        prob = build_nonlinear_interval_problem(
            obs, z0_plus, z0_minus, tspan; x0=z0
        )
    end

    sol = DifferentialEquations.solve(prob, solver)
    return sol
end
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
    # K::Union{Nothing, Vector},
    λ_vals::Tuple{Float64, Float64},
    tspan:: Tuple{Float64, Float64},
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
    tspan::Tuple{Float64, Float64};
    x0::Union{Nothing, Vector} = nothing,
    f_true::Union{Nothing, Vector} = nothing,
    saveat::Union{Nothing, Real, AbstractVector{<:Real}} = nothing,
    solver = Tsit5()
)
    A = sys.A
    C = sys.C
    n = sys.n

    validate_initial_bounds(x0_minus, x0_plus)

    A_minus_KC = A - K * reshape(C, 1, :)
    transformed = false

    if monotone_dynamic(A_minus_KC)
        obs = IntervalObserver(sys, K, f_plus, f_minus)
        prob = build_nonlinear_interval_problem(
            obs, x0_plus, x0_minus, tspan; x0 = x0, f_true = f_true
        )
    else
        transformed = true

        # Change of coordinates z = Mx
        M, M_inv, D = diagonalizing_change_of_basis(A, C, K)
        # @info "Diagonalizing change of basis matrix M: $M"
        # @info "Inverse of M: $M_inv"
        z0_minus, z0_plus, z0 = compute_bounds_in_new_basis(M, x0_minus, x0_plus, x0)
        # @info "= "^20
        # @info "Transformed initial conditions:"
        # @info "z0_minus: $z0_minus"
        # @info "z0_plus: $z0_plus"
        # @info "z0: $z0"
        # @info "= "^20
        A_z = M * A * M_inv
        C_z = vec((C') * M_inv)

        f_plus_z  = transform_function_vector(f_plus,  M)
        f_minus_z = transform_function_vector(f_minus, M)
        f_true_z  = isnothing(f_true) ? nothing : transform_function_vector(f_true, M)

        K_z = M * K

        new_sys = NonLinearSystem(
            A_z, C_z, f_plus_z, f_minus_z;
            check_metzler = false
        )

        obs = IntervalObserver(new_sys, K_z, f_plus_z, f_minus_z)

        prob = build_nonlinear_interval_problem(
            obs, z0_plus, z0_minus, tspan; x0 = z0, f_true = f_true_z
        )
    end

    sol_raw = isnothing(saveat) ?
        DifferentialEquations.solve(prob, solver) :
        DifferentialEquations.solve(prob, solver; saveat = saveat)
    # @info "Solution obtained sol: $sol_raw with status: $(sol_raw.retcode)"
    Z = hcat(sol_raw.u...)
    # @show size(Z)
    # @show Z[1:n, 1]
    # @show Z[1:n, end]
    num_states = size(Z, 1)

    if num_states == 3n
        state_mid = get_state(Z, n)
        lower     = get_lower(Z, n)
        upper     = get_upper(Z, n)
    elseif num_states == 2n
        state_mid = nothing
        upper     = get_upper_nonlinear(Z, n)
        lower     = get_lower_nonlinear(Z, n)
    else
        error("Unexpected state dimension: got $num_states, expected $(2n) or $(3n)")
    end

    if transformed
        xl, xu, x_true = compute_bounds_in_new_basis(M_inv, lower, upper, state_mid)
        @info "Transformed bounds:"

        @show x_true[:, 1]
        @show x_true[:, end]
    else
        xl, xu, x_true = lower, upper, state_mid
    end

    X = isnothing(x_true) ? vcat(xu, xl) : vcat(x_true, xu, xl)
    u = [X[:, k] for k in 1:size(X, 2)]

    return (t = sol_raw.t, u = u)
end

function solve(
    sys::NonLinearSystem,
    λ_vals::Tuple{Float64, Float64},
    f_plus::Vector,
    f_minus::Vector,
    x0_plus::Vector,
    x0_minus::Vector,
    tspan::Tuple{Float64, Float64};
    x0::Union{Nothing, Vector} = nothing,
    f_true::Union{Nothing, Vector} = nothing,
    δ::Float64 = 0.5,
    saveat::Union{Nothing, AbstractVector{<:Real}} = nothing,
    num_saveat::Int = 1001,
    solver = Tsit5()
)
    n = sys.n
    desired_poles = generate_poles(λ_vals, n; δ = δ)

    @assert num_saveat ≥ 2 "num_saveat must be ≥ 2"
    common_t = isnothing(saveat) ?
        collect(range(tspan[1], tspan[2], length = num_saveat)) :
        collect(Float64.(saveat))

    results = Vector{IntervalObserverSolution}(undef, length(desired_poles))

    Threads.@threads for k in eachindex(desired_poles)
        poles = desired_poles[k]

        K = positive_interval_gain(sys, desired_poles = poles)

        sol_observer = solve(
            sys,
            K,
            f_plus,
            f_minus,
            x0_plus,
            x0_minus,
            tspan;
            x0 = x0,
            f_true = f_true,
            saveat = common_t,
            solver = solver,
        )

        results[k] = IntervalObserverSolution(
            sol_observer.t,
            sol_observer.u;
            label = "Desired poles: $(poles)",
            show_true = (k == 1),
        )
    end

    return intersect_solutions(results)
end
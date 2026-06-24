using DifferentialEquations
using Base.Threads

function solve(
    sys::LinearSystem,
    x0::Vector,
    xl0::Vector,
    xu0::Vector,
    K::Vector,
    tspan::Tuple{Float64, Float64},
    solver = Tsit5(),
)
    return solve(
        sys,
        K,
        x0,
        xl0,
        xu0,
        tspan,
        solver,
    )
end

# --------------------------------------------------
# 1. Core nonlinear solver: user provides K
# --------------------------------------------------

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
    solver = Tsit5(),
    M::Union{Nothing, AbstractMatrix} = nothing,
    M_inv::Union{Nothing, AbstractMatrix} = nothing,
    label::String = "single observer",
    show_true::Bool = true,
    )
    A = sys.A
    C = sys.C
    n = sys.n

    validate_initial_bounds(x0_minus, x0_plus)

    A_minus_KC = A - K * reshape(C, 1, :)
    
    if isnothing(M)
        transformed = !monotone_dynamic(A_minus_KC)

        if transformed
            M, M_inv, _ = diagonalizing_change_of_basis(A, C, K)
        end
    else
        transformed = true
        @assert !isnothing(M_inv) "If M is provided, M_inv must also be provided."
    end

    if !transformed
        obs = IntervalObserver(sys, K)

        prob = build_nonlinear_interval_problem(
            obs,
            x0_plus,
            x0_minus,
            tspan;
            x0 = x0,
            f_true = f_true,
        )
    else
        z0_minus, z0_plus, z0 =
            compute_bounds_in_new_basis(M, x0_minus, x0_plus, x0)

        A_z = M * A * M_inv
        C_z = vec((C') * M_inv)
        K_z = M * K

        f_plus_z  = transform_function_vector(f_plus, M)
        f_minus_z = transform_function_vector(f_minus, M)
        f_true_z  = isnothing(f_true) ? nothing : transform_function_vector(f_true, M)

        new_sys = NonLinearSystem(
            A_z,
            C_z,
            f_plus_z,
            f_minus_z;
            check_metzler = false,
        )

        obs = IntervalObserver(new_sys, K_z)

        prob = build_nonlinear_interval_problem(
            obs,
            z0_plus,
            z0_minus,
            tspan;
            x0 = z0,
            f_true = f_true_z,
        )
    end

    sol_raw = _solve_ode(prob, solver, saveat)

    lower, upper, x_true = _extract_solution_blocks(sol_raw, n)

    if transformed
        lower, upper, x_true =
            compute_bounds_in_new_basis(M_inv, lower, upper, x_true)
    end

    return _build_solution(
        sol_raw.t,
        lower,
        upper,
        x_true;
        label = label,
        show_true = show_true,
    )
end

# --------------------------------------------------
# 2. Scalar σ: companion/Vandermonde method
# --------------------------------------------------

function _solve_with_sigma(
    sys::NonLinearSystem,
    σ::Real,
    f_plus::Vector,
    f_minus::Vector,
    x0_plus::Vector,
    x0_minus::Vector,
    tspan::Tuple{Float64, Float64};
    x0::Union{Nothing, Vector} = nothing,
    f_true::Union{Nothing, Vector} = nothing,
    saveat::Union{Nothing, Real, AbstractVector{<:Real}} = nothing,
    solver = Tsit5(),
    )
    M, M_inv, desired_poles = sigma_change_of_basis(sys, σ)

    K = positive_interval_gain(
        sys;
        desired_poles = desired_poles,
    )

    return solve(
        sys,
        K,
        f_plus,
        f_minus,
        x0_plus,
        x0_minus,
        tspan;
        x0 = x0,
        f_true = f_true,
        saveat = saveat,
        solver = solver,
        M = M,
        M_inv = M_inv,
        label = "σ = $(σ), poles = $(desired_poles)",
        show_true = true,
    )
end

# --------------------------------------------------
# 3. One-value tuple: redirect to scalar σ
# --------------------------------------------------

function solve(
    sys::NonLinearSystem,
    λ_vals::Tuple{<:Real},
    f_plus::Vector,
    f_minus::Vector,
    x0_plus::Vector,
    x0_minus::Vector,
    tspan::Tuple{Float64, Float64};
    x0::Union{Nothing, Vector} = nothing,
    f_true::Union{Nothing, Vector} = nothing,
    saveat::Union{Nothing, Real, AbstractVector{<:Real}} = nothing,
    solver = Tsit5(),
    )
    return solve(
        sys,
        λ_vals[1],
        f_plus,
        f_minus,
        x0_plus,
        x0_minus,
        tspan;
        x0 = x0,
        f_true = f_true,
        saveat = saveat,
        solver = solver,
    )
end

# --------------------------------------------------
# 4. Interval λ_vals: family of observers + intersection
# --------------------------------------------------

function solve(
    sys::NonLinearSystem,
    λ_vals::Tuple{<:Real, <:Real},
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
    solver = Tsit5(),
    )
    n = sys.n
    desired_poles = generate_poles(λ_vals, n; δ = δ)
    common_t = _common_time_grid(tspan, saveat, num_saveat)

    results = Vector{IntervalObserverSolution}(undef, length(desired_poles))

    @threads for k in eachindex(desired_poles)
        poles = desired_poles[k]

        K = positive_interval_gain(
            sys;
            desired_poles = poles,
        )

        results[k] = solve(
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
            label = "Desired poles: $(poles)",
            show_true = (k == 1),
        )
    end

    return intersect_solutions(results)
end


# --------------------------------------------------
# Scalar σ solve:
# multi-σ parallel intersection
# --------------------------------------------------

function solve(
    sys::NonLinearSystem,
    σ::Real,
    f_plus::Vector,
    f_minus::Vector,
    x0_plus::Vector,
    x0_minus::Vector,
    tspan::Tuple{Float64, Float64};
    x0::Union{Nothing, Vector} = nothing,
    f_true::Union{Nothing, Vector} = nothing,
    σ_min::Float64 = 0.1,
    num_sigma::Int = 5,
    saveat::Union{Nothing, Real, AbstractVector{<:Real}} = nothing,
    num_saveat::Int = 1001,
    solver = Tsit5(),
)
    σ_values = generate_sigma_family(
        σ;
        σ_min = σ_min,
        num = num_sigma,
    )

    common_t = _common_time_grid(
        tspan,
        saveat,
        num_saveat,
    )

    results = Vector{IntervalObserverSolution}(
        undef,
        length(σ_values),
    )

    @threads for k in eachindex(σ_values)

        σk = σ_values[k]

        results[k] = _solve_with_sigma(
            sys,
            σk,
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
    end

    return intersect_solutions(results)
end
using DifferentialEquations
using Base.Threads

"""
    solve(sys::LinearSystem, x0, xl0, xu0, K, tspan, solver=Tsit5())

Solve a linear interval observer problem.

Solves the coupled ODE system:
  ẋ = A*x
  ẋ⁺ = A*x⁺ + K*(y - C*x⁻)
  ẋ⁻ = A*x⁻ + K*(y - C*x⁺)

where y = C*x is the measurement.

# Arguments
- `sys::LinearSystem`: The linear system
- `x0::Vector`: Initial state (dimension n)
- `xl0::Vector`: Initial lower bound estimate (dimension n)
- `xu0::Vector`: Initial upper bound estimate (dimension n)
- `K::Vector`: Observer gain (dimension n)
- `tspan::Tuple{Float64, Float64}`: Time span (t0, tf)
- `solver`: ODE solver algorithm (default: Tsit5())

# Returns
- `IntervalObserverSolution`: Solution object containing time history and state estimates

# Example
```julia
sys = LinearSystem(A, C)
x0 = [1.0, 0.0]
xl0 = [0.8, -0.2]
xu0 = [1.2, 0.2]
K = positive_interval_gain(sys)
sol = solve(sys, x0, xl0, xu0, K, (0.0, 10.0))
plot(sol)
```
"""
function solve(
    sys::LinearSystem,
    x0::Vector,
    xl0::Vector,
    xu0::Vector,
    K::Vector,
    tspan::Tuple{Float64, Float64},
    solver = Tsit5(),
)
    n = sys.n
    
    # Initial condition state vector: [x; xu; xl]
    # This order matches the Linear_syst_int_obs_ode! function's view structure
    z0 = vcat(x0, xu0, xl0)
    
    # ODE problem parameters
    p = (sys.A, sys.C, K, n)
    
    # Build and solve ODE problem
    prob = ODEProblem(Linear_syst_int_obs_ode!, z0, tspan, p)
    sol_raw = _solve_ode(prob, solver, nothing)
    
    # Extract solution blocks: the ODE returns [x; xu; xl]
    # _extract_solution_blocks expects that layout
    lower, upper, x_true = _extract_solution_blocks(sol_raw, n)
    
    # Build solution object: reconstructs as [x; xu; xl]
    return _build_solution(
        sol_raw.t,
        lower,
        upper,
        x_true;
        label = "linear observer",
        show_true = true,
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

"""
    _solve_with_sigma(sys::NonLinearSystem, σ::Real, f_plus::Vector, f_minus::Vector, 
                      x0_plus::Vector, x0_minus::Vector, tspan; ...)

Internal function to solve with a single sigma parameter.

Uses the sigma-based companion form change of basis to design the observer gain.

# Arguments
- `sys::NonLinearSystem`: The nonlinear system
- `σ::Real`: Sigma parameter > 1 (controls pole placement scaling)
- `f_plus, f_minus, x0_plus, x0_minus, tspan`: See `solve` documentation
- `x0, f_true, saveat, solver`: Optional parameters

# Returns
- `IntervalObserverSolution`: Solution for this sigma value
"""
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

"""
    solve(sys::NonLinearSystem, λ_vals::Tuple{<:Real}, f_plus::Vector, f_minus::Vector, 
          x0_plus::Vector, x0_minus::Vector, tspan; ...)

Solve nonlinear interval observer with single eigenvalue (one-value tuple).

Redirects to scalar sigma-based solver.

# Arguments
- `λ_vals::Tuple{<:Real}`: Single-element tuple containing eigenvalue
- Other arguments: See main `solve` documentation

# Returns
- `IntervalObserverSolution`: Solution object
"""
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

"""
    solve(sys::NonLinearSystem, λ_vals::Tuple{<:Real, <:Real}, f_plus::Vector, f_minus::Vector, 
          x0_plus::Vector, x0_minus::Vector, tspan; δ=0.5, ...)

Solve nonlinear interval observer with eigenvalue interval (two-value tuple).

Generates a family of observer gains with poles distributed in the interval [λ_min, λ_max],
solves each observer, and returns the intersection of all solutions for improved bounds.

# Arguments
- `λ_vals::Tuple{<:Real, <:Real}`: (λ_min, λ_max) eigenvalue bounds
- `δ::Float64`: Geometric distribution parameter (default: 0.5)
- `num_saveat::Int`: Number of time points (default: 1001)
- Other arguments: See main `solve` documentation

# Returns
- `IntervalObserverSolution`: Intersection of all observer solutions
"""
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


"""
    solve(sys::NonLinearSystem, σ::Real, f_plus::Vector, f_minus::Vector, 
          x0_plus::Vector, x0_minus::Vector, tspan; σ_min=1.1, num_sigma=5, ...)

Solve nonlinear interval observer with multi-sigma parallel intersection.

Generates a family of sigma values from σ_min to σ, solves multiple observers
in parallel, and returns the intersection of all solutions for refined bounds.

# Arguments
- `σ::Real`: Maximum sigma parameter (must be > 1)
- `σ_min::Float64`: Minimum sigma parameter (default: 1.1)
- `num_sigma::Int`: Number of sigma values to try (default: 5)
- `num_saveat::Int`: Number of output time points (default: 1001)
- Other arguments: See main `solve` documentation

# Returns
- `IntervalObserverSolution`: Intersection of all sigma-based observer solutions

# Note
Solutions are computed in parallel using multiple threads.
"""
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
    σ_min::Float64 = 1.1,
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
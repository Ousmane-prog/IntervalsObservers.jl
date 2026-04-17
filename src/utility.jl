"""
    get_state(sol, n::Int)

Extract the true state trajectory `x(t)` from the interval observer solution.

# Arguments
- `sol`: ODE solution from interval observer simulation
- `n::Int`: System dimension

# Returns
- `Matrix`: State trajectory (n×T) where T is the number of time points
"""
function get_state(sol, n::Int)
    return sol[1:n, :]
end

"""
    get_lower(sol, n::Int)

Extract the lower bound trajectory `x⁻(t)` from the interval observer solution.

# Arguments
- `sol`: ODE solution from interval observer simulation
- `n::Int`: System dimension

# Returns
- `Matrix`: Lower bound trajectory (n×T) where T is the number of time points
"""
function get_lower(sol, n::Int)
    return sol[2n+1:3n, :] 
end

"""
    get_upper(sol, n::Int)

Extract the upper bound trajectory `x⁺(t)` from the interval observer solution.

# Arguments
- `sol`: ODE solution from interval observer simulation
- `n::Int`: System dimension

# Returns
- `Matrix`: Upper bound trajectory (n×T) where T is the number of time points
"""
function get_upper(sol, n::Int)
    return sol[n+1:2n, :]
end

"""
    get_upper_nonlinear(sol, n::Int)

Extract the upper bound trajectory `x⁺(t)` from a nonlinear interval observer solution.

# Arguments
- `sol`: ODE solution from nonlinear interval observer simulation
- `n::Int`: System dimension

# Returns
- `Matrix`: Upper bound trajectory (n×T) where T is the number of time points
"""
function get_upper_nonlinear(sol, n::Int)
    return sol[1:n, :]
end

"""
    get_lower_nonlinear(sol, n::Int)

Extract the lower bound trajectory `x⁻(t)` from a nonlinear interval observer solution.

# Arguments
- `sol`: ODE solution from nonlinear interval observer simulation
- `n::Int`: System dimension

# Returns
- `Matrix`: Lower bound trajectory (n×T) where T is the number of time points
"""
function get_lower_nonlinear(sol, n::Int)
    return sol[n+1:2n, :]
end

function transform_function_vector(f_vec::Vector, M::AbstractMatrix)
    n = length(f_vec)

    return [
        (t, y) -> sum(M[i, j] * f_vec[j](t, y) for j in 1:n)
        for i in 1:n
    ]
end
# ============================================================================
# Observer Gain Constants
# ============================================================================

"""
Default gain value for states with positive measurement."""
const DEFAULT_GAIN_VALUE = 0.1

"""
Scaling factor for desired poles placement."""
const POLE_PLACEMENT_SCALE = 0.5

# ============================================================================
# Helper Functions for Observer Gain Computation
# ============================================================================

"""
    _find_non_monotone_entries(M::Matrix{<:Real}) -> Vector{Tuple{Int, Int, Float64}}

Identify all non-positive off-diagonal entries in matrix M that violate monotonicity.

# Returns
- `Vector{Tuple{Int, Int, Float64}}`: Tuples of (row, col, value) for non-positive entries
"""
function _find_non_monotone_entries(M::Matrix{<:Real})
    non_positive = Tuple{Int, Int, Float64}[]
    n = size(M, 1)
    
    @inbounds for i in 1:n
        for j in 1:n
            if i != j && M[i, j] < 0
                push!(non_positive, (i, j, M[i, j]))
            end
        end
    end
    
    return non_positive
end

"""
    _is_monotone_dynamic(M::Matrix{<:Real}) -> Bool

Check if matrix M has a monotone (Metzler) structure: all off-diagonal entries are non-negative.

# Arguments
- `M::Matrix{<:Real}`: Matrix to check

# Returns
- `Bool`: true if matrix is monotone, false otherwise
"""
function _is_monotone_dynamic(M::Matrix{<:Real})
    return isempty(_find_non_monotone_entries(M))
end

"""
    _compute_default_gain(C::Vector{<:Real}, n::Int) -> Vector{<:Real}

Compute a simple diagonal gain based on measurement availability.

Sets K[i] = DEFAULT_GAIN_VALUE if C[i] > 0 (state i is measured), 0 otherwise.

# Arguments
- `C::Vector{<:Real}`: Measurement vector
- `n::Int`: System dimension

# Returns
- `Vector{<:Real}`: Gain vector K
"""
function _compute_default_gain(C::Vector{<:Real}, n::Int)
    K = zeros(length(C))
    
    @inbounds for i in 1:n
        if C[i] >= 0
            K[i] = DEFAULT_GAIN_VALUE
        end
    end
    
    if all(k == 0 for k in K)
        throw(UnobservableMeasurementError(
            "No states are measured: C vector has no positive entries. " *
            "At least one C[i] > 0 is required for interval observer design.",
            C
        ))
    end
    
    return K
end

"""
    _validate_desired_poles(desired_poles::Vector{<:Real}, n::Int) -> Bool

Validate that desired poles are distinct and correct in number.

# Arguments
- `desired_poles::Vector{<:Real}`: Desired eigenvalues for error dynamics
- `n::Int`: System dimension (must equal length of desired_poles)

# Throws
- `InvalidDesiredPolesError`: If poles are not distinct or have wrong dimension

# Returns
- `Bool`: true if poles are valid
"""
function _validate_desired_poles(desired_poles::Vector{<:Real}, n::Int)
    # Check dimension
    if length(desired_poles) != n
        throw(InvalidDesiredPolesError(
            "Desired poles dimension ($(length(desired_poles))) must match system dimension ($n)",
            desired_poles,
            "dimension mismatch"
        ))
    end
    
    # Check that poles are distinct (no repeated values)
    unique_poles = length(unique(desired_poles))
    if unique_poles != n
        throw(InvalidDesiredPolesError(
            "Desired poles must be distinct. Found $(n - unique_poles) repeated value(s) in poles: $(desired_poles)",
            desired_poles,
            "repeated poles"
        ))
    end
    
    return true
end

"""
    _compute_placed_gain(A::Matrix{<:Real}, C::Vector{<:Real}, desired_poles::Vector{<:Real}, n::Int) -> Vector{<:Real}

Compute observer gain using pole placement for error dynamics A - K*C.

This function uses pole placement (Ackermann's formula) to compute observer gain K such that
the eigenvalues (poles) of the error dynamics matrix (A - K*C) are placed at the desired locations.

When the user specifies `desired_poles = [λ₁, λ₂, ..., λₙ]`, this function finds K such that:
  eig(A - K*C) = [λ₁, λ₂, ..., λₙ]

Negative eigenvalues ensure exponential convergence of the observation error.

# Arguments
- `A::Matrix{<:Real}`: System matrix
- `C::Vector{<:Real}`: Measurement vector
- `desired_poles::Vector{<:Real}`: Desired eigenvalues (poles) for the error dynamics A - K*C
  These are the desired values of eig(A - K*C). Must be distinct real numbers.
- `n::Int`: System dimension

# Returns
- `Vector{<:Real}`: Gain vector K

# Throws
- `InvalidDesiredPolesError`: If poles are not distinct or have wrong dimension
- `NonMonotoneDynamicsError`: If the resulting A - K*C is not monotone (Metzler)
"""
function _compute_placed_gain(A::Matrix{<:Real}, C::Vector{<:Real}, desired_poles::Vector{<:Real}, n::Int)
    # Validate poles first
    _validate_desired_poles(desired_poles, n)
    
    # The desired_poles are the eigenvalues we want for the error dynamics matrix (A - K*C)
    # scaled_poles = POLE_PLACEMENT_SCALE .* desired_poles
    L = place(A', reshape(C, :, 1), desired_poles)
    K = vec(L)
    
    # Check monotonicity of error dynamics matrix (A - K*C)
    # Note: The eigenvalues of (A - K*C) should equal the desired_poles we specified
    A_minus_KC = A - K * reshape(C, 1, :)
    
    # if !_is_monotone_dynamic(A_minus_KC)
    #     # non_monotone = _find_non_monotone_entries(A_minus_KC)
    #     T, M = diagonalize_matrix(A_minus_KC)
    # end
    
    return K
end

"""
    diagonalizing_change_of_basis(A, C, K)

Compute the change of basis matrix M such that M(A - K*C)M⁻¹ is diagonal.

retruns
   M         : the change of basis matrix
   M_inv     : the inverse of the change of basis matrix
   D         : the diagonal matrix of eigenvalues
"""
function diagonalizing_change_of_basis(A, C, K)
    A_minus_KC = A - K * reshape(C, 1, :)
    F = eigen(A_minus_KC)
    V = F.vectors
    D = Diagonal(F.values)
    M_inv = V 
    M = inv(V)
    return M, M_inv, D    
end


function transform_system(A, C, f, M)
    M_inv = inv(M)
    A_new = M * A * M_inv
    C_new = vec(C * M_inv)
end

function diagonalize_matrix(M)
    F = eigen(M)
    return F.vectors, Diagonal(F.values)
end


# function transform_interval_bounds(
#     P::AbstractMatrix,
#     lower::Union{AbstractMatrix, AbstractVector},
#     upper::Union{AbstractMatrix, AbstractVector},
#     state::Union{AbstractMatrix, AbstractVector, Nothing} = nothing
# )
#     y_lower_raw = P * lower
#     y_upper_raw = P * upper

#     y_lower = min.(y_lower_raw, y_upper_raw)
#     y_upper = max.(y_lower_raw, y_upper_raw)

#     y_state = isnothing(state) ? nothing : P * state

#     return y_lower, y_upper, y_state
# end

function generate_extreme_points(lower::AbstractVector, upper::AbstractVector)
    function gen(i)
        if i > 1
            E1 = gen(i - 1)
            nb_col1 = size(E1, 2)

            top = vcat(E1, lower[i] .* ones(1, nb_col1))
            bottom = vcat(E1, upper[i] .* ones(1, nb_col1))

            return hcat(top, bottom)
        else
            return reshape([lower[1], upper[1]], 1, 2)
        end
    end

    return gen(length(lower))
end


function compute_bounds_in_new_basis(
    M::AbstractMatrix,
    x_minus::Union{AbstractVector, AbstractMatrix},
    x_plus::Union{AbstractVector, AbstractMatrix},
    x::Union{AbstractVector, AbstractMatrix, Nothing} = nothing
)
    # Single box
    if x_minus isa AbstractVector && x_plus isa AbstractVector
        E = generate_extreme_points(x_minus, x_plus)
        EZ = M * E

        z_minus = vec(minimum(EZ, dims=2))
        z_plus  = vec(maximum(EZ, dims=2))
        z       = isnothing(x) ? nothing : M * x

        return z_minus, z_plus, z
    end

    # Several boxes stored columnwise
    @assert x_minus isa AbstractMatrix && x_plus isa AbstractMatrix
    @assert size(x_minus) == size(x_plus)
    @assert isnothing(x) || size(x) == size(x_minus)

    n_out = size(M, 1)
    n_cols = size(x_minus, 2)

    z_minus = Matrix{eltype(x_minus)}(undef, n_out, n_cols)
    z_plus  = Matrix{eltype(x_plus)}(undef, n_out, n_cols)
    z       = isnothing(x) ? nothing : Matrix{eltype(x)}(undef, n_out, n_cols)

    for k in 1:n_cols
        E = generate_extreme_points(x_minus[:, k], x_plus[:, k])
        EZ = M * E

        z_minus[:, k] = vec(minimum(EZ, dims=2))
        z_plus[:, k]  = vec(maximum(EZ, dims=2))

        if !isnothing(x)
            z[:, k] = M * x[:, k]
        end
    end

    return z_minus, z_plus, z
end
# ============================================================================
# Main Observer Gain Function
# ============================================================================

"""
    positive_interval_gain(sys::Union{LinearSystem, NonLinearSystem}; 
                          desired_poles::Union{Vector{<:Real}, Nothing}=nothing) -> Vector{<:Real}

Compute a positive interval observer gain `K` for the given system.

## Modes of Operation

### Default Mode (no desired_poles)
Uses a simple diagonal gain: K[i] = DEFAULT_GAIN_VALUE if C[i] > 0, else 0.
This is parameter-free and works when C has positive entries.

### Pole Placement Mode (desired_poles provided)
Designs K via pole placement to place the eigenvalues (poles) of error dynamics A - K*C 
at the user-specified locations. The user provides the desired eigenvalues, and this 
function computes K such that: eig(A - K*C) = desired_poles

Automatically scales poles by POLE_PLACEMENT_SCALE (0.5) for numerical stability.
**Important**: Desired poles must be distinct (no repeated values).
Validates that the resulting A - K*C remains monotone (Metzler).

# Arguments
- `sys::Union{LinearSystem, NonLinearSystem}`: The system for which to compute the gain
- `desired_poles::Union{Vector{<:Real}, Nothing}`: 
  - `nothing`: Use default diagonal gain (default)
  - `Vector{<:Real}`: Desired eigenvalues (poles) for the error dynamics matrix A - K*C.
    **These are the characteristic polynomial roots you want:** λ ∈ ℝ such that eig(A - K*C) = [λ₁, λ₂, ..., λₙ]
    Negative values ensure exponential decay of observation error.
    **Must be:** distinct, real-valued, same length as system dimension

# Returns
- `Vector{<:Real}`: Observer gain vector K of dimension n

# Throws
- `UnobservableMeasurementError`: If no states are measured (all C[i] ≤ 0)
- `InvalidDesiredPolesError`: If poles are not distinct or have wrong dimension
- `NonMonotoneDynamicsError`: If pole placement yields non-monotone dynamics

# Example
```julia
sys = LinearSystem(A, C)

# Use default gain
K1 = positive_interval_gain(sys)

# Use pole placement with desired eigenvalues of A - K*C
K2 = positive_interval_gain(sys, desired_poles=[-1.0, -2.0])
```

# Notes
- The "poles" are the eigenvalues of the error dynamics matrix: poles = eig(A - K*C)
- The returned gain K is used in the observer dynamics: dx̂/dt = Ax̂ + K(C(x - x̂))
- For a positive interval observer, all off-diagonal entries of A - K*C must be non-negative (Metzler)
- Pole placement automatically scales desired poles by POLE_PLACEMENT_SCALE for numerical stability
"""
function positive_interval_gain(sys::Union{LinearSystem, NonLinearSystem}; 
                               desired_poles::Union{Vector{<:Real}, Nothing}=nothing)
    n = sys.n
    C = sys.C
    A = sys.A
    
    if desired_poles === nothing
        # Use simple default gain
        K = _compute_default_gain(C, n)
    else
        # Use pole placement
        K = _compute_placed_gain(A, C, desired_poles, n)
    end
    
    return K
end

"""
    monotone_dynamic(M::Matrix{<:Real}) -> Bool

Check if the matrix M is monotone (Metzler): all off-diagonal entries are non-negative.

# Arguments
- `M::Matrix{<:Real}`: Matrix to check

# Returns
- `Bool`: true if matrix is monotone, false otherwise
"""
function monotone_dynamic(M::Matrix{<:Real})
    return _is_monotone_dynamic(M)
end

function create_collection(λ_vals::Tuple{Float64, Float64}, n::Integer)
    λ_min, λ_max = Float64.(λ_vals)
    @assert λ_min < 0 "λ_min must be negative for stability"
    @assert λ_max < 0 "λ_max must be negative for stability"
    @assert λ_min < λ_max "λ_min must be less than λ_max"

    return collect(range(λ_min, λ_max; length=n))
end

function generate_poles_geometric(λ::Float64, n::Integer; δ=0.5)
    @assert λ < 0 "λ must be negative for stability"
    @assert 0 < δ < 1 "δ must be in (0, 1)"

    return [λ * (1 - δ)^k for k in 0:(n - 1)]
end

function generate_poles(λ_vals::Tuple{Float64, Float64}, n::Integer; δ=0.5)
    pole_collection = create_collection(λ_vals, n)
    return [generate_poles_geometric(λ, n; δ = δ) for λ in pole_collection]
end

function solution_to_matrix(sol)
    return reduce(hcat, sol.u)
end

function intersect_solutions(results::Vector{IntervalObserverSolution})
    @assert !isempty(results) "Cannot intersect an empty collection of solutions."

    Z0 = solution_to_matrix(results[1])
    num_states = size(Z0, 1)
    nt = size(Z0, 2)

    has_true = iszero(num_states % 3)

    if has_true
        n = num_states ÷ 3
        x_true = Z0[1:n, :]
        upper_int = copy(Z0[n+1:2n, :])
        lower_int = copy(Z0[2n+1:3n, :])
    else
        n = num_states ÷ 2
        x_true = nothing
        upper_int = copy(Z0[1:n, :])
        lower_int = copy(Z0[n+1:2n, :])
    end

    for k in 2:length(results)
        @assert results[k].t == results[1].t "All solutions must share the same time grid."

        Z = solution_to_matrix(results[k])
        @assert size(Z, 1) == num_states
        @assert size(Z, 2) == nt

        if has_true
            upper_k = Z[n+1:2n, :]
            lower_k = Z[2n+1:3n, :]
        else
            upper_k = Z[1:n, :]
            lower_k = Z[n+1:2n, :]
        end

        upper_int .= min.(upper_int, upper_k)
        lower_int .= max.(lower_int, lower_k)
    end

    if any(lower_int .> upper_int)
        @warn "Empty intersection detected for some states/times."
    end

    X = isnothing(x_true) ? vcat(upper_int, lower_int) : vcat(x_true, upper_int, lower_int)
    u = [X[:, k] for k in 1:size(X, 2)]

    return IntervalObserverSolution(
        results[1].t,
        u,
        label = "intersection",
        show_true = true,
    )
end
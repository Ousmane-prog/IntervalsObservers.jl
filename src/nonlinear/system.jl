"""
    NonLinearSystem{T<:Real, F}

Represents a nonlinear control-affine system for interval observer design.

The system dynamics are:
  ẋ(t) = A*x(t) + f(t, y(t))
  y(t) = C*x(t)

where:
- x(t) ∈ ℝⁿ is the state
- A ∈ ℝⁿˣⁿ is the linear part
- f(t,y) is the nonlinear input that depends on time and measurement y
- C ∈ ℝ¹ˣⁿ is the measurement vector

For positive interval observer design, f must satisfy:
  f┈(t, y) ≤ f(t, y) ≤ f↑(t, y)

where f┈ (f_minus) and f↑ (f_plus) bound the true nonlinearity.

# Fields
- `A::Matrix{T}`: Linear system matrix (n×n)
- `C::Vector{T}`: Measurement vector (length n)
- `n::Int`: System dimension
- `obs::Matrix{T}`: Observability matrix
- `observable::Bool`: Indicates system is observable
- `is_metzler::Bool`: Indicates A is Metzler
- `f_plus::Vector{F}`: Upper bounding functions [f↑₁, ..., f↑ₙ], each signature: (t, y) -> ℝ
- `f_minus::Vector{F}`: Lower bounding functions [f┈₁, ..., f┈ₙ], each signature: (t, y) -> ℝ
- `check_metzler::Bool`: Flag for Metzler matrix checking

# Constructor
```julia
NonLinearSystem(A, C, f_plus, f_minus; check_metzler=true)
```

# Throws
- `DimensionMismatchError`: If A, C, f_plus, or f_minus have inconsistent dimensions
- `NonObservableSystemError`: If system is not observable
- `NonMetzlerMatrixError`: If A is not Metzler (when check_metzler=true)
"""
struct NonLinearSystem{T<:Real, F}
    A::Matrix{T}
    C::Vector{T}
    n::Int
    obs::Matrix{T}
    observable::Bool
    is_metzler::Bool
    f_plus::Vector{F}
    f_minus::Vector{F}
    check_metzler::Bool

    function NonLinearSystem(
        A::Matrix{T},
        C::Vector{T},
        f_plus::Vector{F},
        f_minus::Vector{F};
        check_metzler::Bool = true
    ) where {T<:Real, F}

        n = validate_system_dimensions(A, C)

        length(f_plus) == n || error("f_plus must contain $n functions")
        length(f_minus) == n || error("f_minus must contain $n functions")

        obs = compute_observability_matrix(A, C, n)
        # observable = true
        observable = check_observability_matrix(obs, n)

        is_metzler =
            check_metzler ? (check_Metzler_Matrix(A); true) : false

        return new{T, F}(A, C, n, obs, observable, is_metzler, f_plus, f_minus, check_metzler)
    end
end

"""
    IntervalObserver{T<:Real, F}

Take an observable nonlinear system 'sys' and an observer gain vector 'K' to create an interval observer object.

An interval observer consists of a pair:


  ẋ⁺(t) = A*x(t) + f⁺(t, y(t)) + K*(y(t) - C*x⁻(t))
  ẋ⁻(t) = A*x(t) + f⁻(t, y(t)) + K*(y(t) - C*x⁺(t))

where x⁺ and x⁻ are the upper and lower state estimates respectively.

# Fields
- `sys::NonLinearSystem{T, F}`: The underlying nonlinear system
- `K::Vector{T}`: Observer gain vector (dimension n)

# Example
```julia
sys = NonLinearSystem(A, C, f_plus, f_minus)
K = positive_interval_gain(sys)
obs = IntervalObserver(sys, K)
```
"""
struct IntervalObserver{T<:Real, F}
    sys::NonLinearSystem{T, F}
    K::Vector{T}
end


"""
    desired_polynomial(roots::Vector{T}) where T

Compute polynomial coefficients from its roots.

Given the roots of a polynomial, this function constructs the polynomial
and returns its coefficients in descending order of degree:

  p(s) = (s - r₁)(s - r₂)...(s - rₙ)

# Arguments
- `roots::Vector{T}`: Vector of polynomial roots

# Returns
- `Vector{T}`: Polynomial coefficients in decreasing degree order
"""
function desired_polynomial(roots::Vector{T}) where T

    p = fromroots(roots)
    coeffs_p = coeffs(p)

    return coeffs_p
end


"""
    observable_canonical_form(sys::NonLinearSystem, P)

Transform system matrix A to observable canonical form using the change of basis matrix P.

Computes the transformed system matrix under the change of variables z = P⁻¹*x:
  ẋ = A*x  ⇒  ẋz = (P⁻¹*A*P)*z

# Arguments
- `sys::NonLinearSystem`: The nonlinear system
- `P::Matrix`: Change of basis matrix

# Returns
- `Matrix`: Transformed system matrix in the new basis
"""
function observable_canonical_form(sys::NonLinearSystem, P)

    A = sys.A

    NA = P \ (A * P) 

    return NA
end

using LinearAlgebra

"""
    check_observability_matrix(obs::Matrix, n::Int)

Verify that the observability matrix has full rank.

Throws an error if the observability matrix doesn't have rank n, indicating
the system is not observable.

# Arguments
- `obs::Matrix`: The observability matrix
- `n::Int`: Expected rank (system dimension)

# Returns
- `Bool`: true if observable, throws otherwise

# Throws
- `ErrorException`: If rank(obs) ≠ n
"""
function check_observability_matrix(obs, n)
    rank(obs) == n || error("The pair (A, C) is not observable.")
    return true
end

"""
    interval_observer_gain(sys::NonLinearSystem, roots::Vector)

Compute an interval observer gain via canonical form transformation.

This function:
1. Computes a positive interval gain for the original system
2. Uses this gain as a change of basis matrix to transform A to observable canonical form
3. Designs the gain G in the canonical coordinates using desired eigenvalues (roots)
4. Transforms the gain back to the original coordinates

# Arguments
- `sys::NonLinearSystem`: The nonlinear system
- `roots::Vector`: Desired eigenvalues for the error dynamics

# Returns
- `Vector`: Observer gain vector G in the original state space
"""
function interval_observer_gain(sys::NonLinearSystem, roots::Vector)

    n = sys.n
    A = sys.A

    P = positive_interval_gain(sys)
    NA = observable_canonical_form(sys, P)
   
    rho = coeffs(fromroots(roots))

    # Ensure correct length
    rho = rho[1:n]

    G = P * (-rho .- NA[:, end])

    return G
end
using LinearAlgebra
using ControlSystems

"""
    LinearSystem{T<:Real}

Represents a linear time-invariant system for interval observer design.

A linear system is defined by:
  ẋ(t) = A*x(t)
  y(t) = C*x(t)

where:
- x(t) ∈ ℝⁿ is the system state
- A ∈ ℝⁿˣⁿ is the system matrix
- C ∈ ℝ¹ˣⁿ is the measurement vector (must be observable)

The system is automatically validated to ensure:
- A is a square matrix
- The system is observable (rank of observability matrix = n)
- A is a Metzler matrix (positive off-diagonal entries, critical for positive interval observers)

# Fields
- `A::Matrix{T}`: System matrix (n×n)
- `C::Vector{T}`: Measurement vector (length n)
- `n::Int`: System dimension
- `observable::Bool`: Indicator that the system is observable
- `is_metzler::Bool`: Indicator that A is Metzler (all off-diagonals are non-negative)

# Constructor
```julia
LinearSystem(A::Matrix{<:Real}, C::Vector{<:Real})
```

# Throws
- `DimensionMismatchError`: If A is not square or C has incorrect length
- `NonObservableSystemError`: If the system is not observable
- `NonMetzlerMatrixError`: If A is not Metzler

# Example
```julia
A = [0.0 1.0; -1.0 -0.5]
C = [1.0, 0.0]
sys = LinearSystem(A, C)
```
"""
struct LinearSystem{T<:Real}
    A::Matrix{T}
    C::Vector{T}
    n::Int
    observable::Bool
    is_metzler::Bool

    function LinearSystem(A::Matrix{T}, C::Vector{T}) where T<:Real
        # n, p = size(A)
        # n == p || error("A must be square matrix, got $(size(A))")
        n = validate_system_dimensions(A, C)
        
        # Check observability criterion: rank([C; CA; CA²; ...; CA^(n-1)]) = n
        M = compute_observability_matrix(A, C, n)
        check_observability(M, n)
        
        observable = true
        check_Metzler_Matrix(A)
        is_metzler = true

        new{T}(A, C, n, observable, is_metzler)
    end 
end


"""
    compute_observability_matrix(A::Matrix{T}, C::Vector{T}, n::Int) where T<:Real

Compute the observability matrix for a linear system.

For a linear system ẋ = Ax, y = Cx, the observability matrix is constructed as:

```
M = [C; C*A; C*A²; ...; C*A^(n-1)]
```

The system is observable if and only if rank(M) = n.

# Arguments
- `A::Matrix{T}`: System matrix (n×n)
- `C::Vector{T}`: Measurement vector (1×n)
- `n::Int`: System dimension

# Returns
- `Matrix{T}`: Observability matrix of size n×n

# Example
```julia
A = [0.0 1.0; -1.0 -0.5]
C = [1.0, 0.0]
n = 2
M = compute_observability_matrix(A, C, n)
rank(M) == n  # Check observability
```
"""
function compute_observability_matrix(A::Matrix{T}, C::Vector{T}, n::Int) where T<:Real
    C_row = reshape(C, 1, n)
    
    M = zeros(T, n, n)
    
    # First row: C
    M[1, :] = C_row
    
    # Subsequent rows: 
    CA = C_row
    @inbounds for i in 2:n
        CA = CA * A
        M[i, :] = CA
        # @info M[i, :]
    end
    
    return M
end



# function positive_interval_gain(sys::Union{LinearSystem, NonLinearSystem})
#     # sys.positive || error("interval observers require a positive system")

#     n = sys.n
#     C = sys.C

#     K = zeros(n)
#     for i in 1:n
#         if C[i] > 0
#             # K[i] = 1.0
#             K[i] = 0
#         end 
#     end
#     return K
# end

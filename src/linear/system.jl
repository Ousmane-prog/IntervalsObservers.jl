using LinearAlgebra
using ControlSystems


struct LinearSystem{T<:Real}
    A::Matrix{T}
    C::Vector{T}
    n::Int
    observable::Bool
    positive::Bool

    function LinearSystem(A::Matrix{T}, C::Vector{T}) where T<:Real
        n, p = size(A)
        n == p || error("A must be square matrix, got $(size(A))")
        
        length(C) == n || error("C length must match A dimensions, got A: $(size(A)), C length: $(length(C))")
        
        # Check observability criterion: rank([C; CA; CA²; ...; CA^(n-1)]) = n
        M = compute_observability_matrix(A, C, n)
        observable = rank(M) == n
        
        # Check if system is positive 
        is_metzler = all(A[i, j] >= 0 || i == j for i in 1:n, j in 1:n)
        C_positive = all(C .>= 0)

        positive = is_metzler && C_positive
        new{T}(A, C, n, observable, positive)
    end 
end


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
    end
    
    return M
end



function positive_interval_gain(sys::LinearSystem)
    sys.positive || error("interval observers require a positive system")

    n = sys.n
    C = sys.C

    K = zeros(n)
    for i in 1:n
        if C[i] > 0
            K[i] = 1.0
        end 
    end
    return K
end



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
    return sol[n+1:2n, :]
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
    return sol[2n+1:3n, :]
end
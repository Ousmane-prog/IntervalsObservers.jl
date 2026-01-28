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
        
        # Check if system is positive (Metzler A and non-negative C)
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
    for i in 2:n
        CA = CA * A
        M[i, :] = CA
    end
    
    return M
end


is_observable(sys::LinearSystem) = sys.observable


function observer_gain(sys::LinearSystem; Q = I)
    sys.observable || error("System is not observable — observer gain does not exist")
    
    A = sys.A
    C = reshape(sys.C, 1, sys.n)

    # Solve dual Riccati equation
    P = care(A', C', Q, 1.0)
    K = -(P*C')
    
    # Return as vector
    return vec(K)
end


function state_observer(sys::LinearSystem; Q = I(sys.n))
    K = observer_gain(sys; Q=Q)

    C_row = reshape(sys.C, 1, sys.n)
    A_obs = sys.A - K * C_row  # Note: (A - KC) not (A + KC)
    B_obs = K

    return A_obs, B_obs, K
end

function upper_observer(sys::LinearSystem; Q = I)
    sys.observable || error("System is not observable — observer does not exist")
    sys.positive || @warn "System is not positive - upper observer may not guarantee bounds"
    
    # Compute gain using CARE
    K = observer_gain(sys; Q=Q)
    
    # Upper observer matrices
    C_row = reshape(sys.C, 1, sys.n)
    A_upper = sys.A - K * C_row
    B_upper = K
    
    # Check if A_upper is Metzler (required for monotone system)
    n = sys.n
    is_metzler = all(A_upper[i, j] >= -1e-10 || i == j for i in 1:n, j in 1:n)
    if !is_metzler
        @warn "Upper observer matrix (A - KC) is not Metzler - bounds may not be guaranteed"
        println("A_upper = ")
        display(A_upper)
    end
    
    return A_upper, B_upper, K
end


function lower_observer(sys::LinearSystem; Q = I)
    sys.observable || error("System is not observable — observer does not exist")
    sys.positive || @warn "System is not positive - lower observer may not guarantee bounds"
    
    # Compute gain using CARE
    K = observer_gain(sys; Q=Q)
    
    # Lower observer matrices (same structure as upper for symmetric design)
    C_row = reshape(sys.C, 1, sys.n)
    A_lower = sys.A - K * C_row
    B_lower = K
    
    # Check if A_lower is Metzler
    n = sys.n
    is_metzler = all(A_lower[i, j] >= -1e-10 || i == j for i in 1:n, j in 1:n)
    if !is_metzler
        @warn "Lower observer matrix (A - KC) is not Metzler - bounds may not be guaranteed"
        println("A_lower = ")
        display(A_lower)
    end
    
    return A_lower, B_lower, K
end


function check_metzler(A::Matrix{T}; tol=1e-10) where T<:Real
    n = size(A, 1)
    return all(A[i, j] >= -tol || i == j for i in 1:n, j in 1:n)
end

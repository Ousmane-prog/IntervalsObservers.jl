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
function get_lower(sol, n::Int)
    return sol[2n+1:3n, :]
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

"""
    positive_interval_gain(sys::Union{LinearSystem, NonLinearSystem}; desired_poles::Union{Vector{<:Real}, Nothing}=nothing)
Compute a positive interval observer gain `K` for the given system. If `desired_poles` are provided, the gain is designed to 
place the eigenvalues of `A-KC` at those locations.
# Arguments
- `sys`: The system (linear or nonlinear) for which to compute the gain
- `desired_poles`: Optional vector of desired eigenvalues for the observer error dynamics
"""
function positive_interval_gain(sys::Union{LinearSystem, NonLinearSystem}; desired_poles::Union{Vector{<:Real}, Nothing}=nothing)

    n = sys.n
    C = sys.C
    A = sys.A
    if desired_poles === nothing
        K = zeros(n)
        for i in 1:n
            if C[i] > 0
                K[i] = 0.1
                # K[i] = 0
            end 
        end
    else
        L = place(A', reshape(C, :, 1), 0.5*desired_poles)
        K = vec(L)
        # 
        A_moinsKC = A - K * reshape(C, 1, :)
        eigvals_moinsKC = eigvals(A_moinsKC)
        println("A-KC: ", A_moinsKC)
        println("Eigenvalues of A - K*C: ", eigvals_moinsKC)
    end
    return K
end

#  ecrire des tests pour les fonctions de utility.jl
# using Test
# using IntervalsObservers



struct NonLinearSystem{T<:Real, F}
    A::Matrix{T}
    C::Vector{T}
    n::Int
    obs::Matrix{T}
    observable::Bool
    is_metzler::Bool
    f_plus::Vector{F}
    f_minus::Vector{F}

    function NonLinearSystem(
        A::Matrix{T},
        C::Vector{T},
        f_plus::Vector{F},
        f_minus::Vector{F}
    ) where {T<:Real, F}

        n = validate_system_dimensions(A, C)

        # Check correct length of nonlinear terms
        # length(f) == n || error("f must contain $n functions")
        length(f_plus) == n || error("f_plus must contain $n functions")
        length(f_minus) == n || error("f_minus must contain $n functions")

        # Observability
        obs = compute_observability_matrix(A, C, n)
        observable = true

        check_Metzler_Matrix(A)
        is_metzler = true

        new{T, F}(A, C, n, obs, observable, is_metzler, f_plus, f_minus)
    end
end

struct IntervalObserver{T<:Real, Fp, Fm}
    sys::NonLinearSystem{T}
    K::Vector{T}
    f_plus::Fp
    f_minus::Fm
end


function desired_polynomial(roots::Vector{T}) where T

    p = fromroots(roots)
    coeffs_p = coeffs(p)

    return coeffs_p
end


function observable_canonical_form(sys::NonLinearSystem, P)

    A = sys.A

    NA = P \ (A * P) 

    return NA
end


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
using LinearAlgebra
using ControlSystems


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

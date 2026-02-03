function Linear_syst_int_obs_ode!(dz, z, p, t)
    A, C, K, n = p

    @inbounds begin

        x =view(z, 1:n)
        xu = view(z, n+1:2n)
        xl = view(z, 2n+1:3n)

        dx = view(dz, 1:n)
        dxu = view(dz, n+1:2n)
        dxl = view(dz, 2n+1:3n)
    end   

    y = dot(C, x)
    ### this part allocates more or less
    dx .= A*x
    dxl .= A*xl + K*(y-dot(C, xu))
    dxu .= A*xu + K*(y-dot(C, xl))
end



# # Innovation terms
#     innov_l = y - Cxu  # y - Cᵀx⁺
#     innov_u = y - Cxl  # y - Cᵀx⁻
    
#     # Compute derivatives using mul! for in-place matrix-vector multiplication
#     mul!(dx, A, x)      # dx = Ax
#     mul!(dxl, A, xl)    # dxl = Axl (will add innovation next)
#     mul!(dxu, A, xu)    # dxu = Axu (will add innovation next)
    
#     # Add innovation terms: dxl += K*innov_l, dxu += K*innov_u
#     @inbounds @simd for i in 1:n
#         dxl[i] += K[i] * innov_l
#         dxu[i] += K[i] * innov_u
#     end
    
#     return nothing
# end
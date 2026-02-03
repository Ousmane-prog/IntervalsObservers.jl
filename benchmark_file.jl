using BenchmarkTools
using Profile
using LinearAlgebra


A = [2.0 3.0 4.0;
    5.0 6.0 7.0;
    2.0 3.0 5.0]
n, p = size(A)
# B = reshape([1 2 3], 3, 1)
B = [1.0, 2.0, 3.0]
C = A*B

function compute_observability_matrix_bis(
    A::Matrix{T},
    C::Vector{T}) where T<:Real

    n, p = size(A)
    C_row = reshape(C,1, n)
    M = zeros(T, n, n)
    CA = C_row
    M[1,:] = CA
    for i in 2:n
        CA = CA*A
        M[i, :] = CA
       
    end
    return M
    
end
# compute_observability_matrix_bis(A, B)
@benchmark display(compute_observability_matrix_bis($A, $B))



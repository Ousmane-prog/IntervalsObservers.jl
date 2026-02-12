struct NonLinearSystem{T<:Real} 
    A::Matrix{T}
    C::Vector{T}
    n::Int 
    f::Function

    function NonLinearSystem(A::Matrix{T}, C::Vector{T}, n::Int, f::Function) where T<:Real
         n = validate_system_dimensions(A, C)
        new{T}(A, C, n, f)
    end

end

macro def(expr)
   quote
        let 
            A, C, f, n 
            $(esc(expr))
            NonLinearSystem(
                A = A,
                C = C,
                n = n,
                f = f)
        end
   end
end
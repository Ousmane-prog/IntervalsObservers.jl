"""
    IntervalObserverError

Abstract type for all interval observer-related errors.
"""

abstract type IntervalObserversError <: Exception end


"""
    NonMetzlerMatrixError

Thrown when attemting to create an interval observer with non metzler Matrix
# Fields 
- 'msg::String': Error message
- 'A::Matrix': The system matrix
"""
struct NonMetzlerMatrixError <: IntervalObserversError
    msg::String
    A::Union{Nothing, Matrix}
end

"""
    check_Metzler_Matrix(A::Matrix{T}) where T

Check if the Matrix A is positive Metzler A ie A[i, j] > 0 if i != j.
Throws NonMetzlerMatrixError with detailed diagnostics if not.
"""
function check_Metzler_Matrix(A::Matrix{T}) where T
    n = size(A, 1)
    
    # Check Metzler
    @inbounds for j in 1:n, i in 1:n
        if i != j && A[i, j] < 0
            throw(NonMetzlerMatrixError(
                "The input Matrix is not Metzler: A[$i,$j] = $(A[i,j])",
                A
            ))
        end
    end

    # # Check C non-negative
    # C_positive = true
    # @inbounds for i in 1:n
    #     if C[i] < 0
    #         C_positive = false
    #         break
    #     end
    # end
    
    # if !is_metzler || !C_positive
    #     throw(NonPositiveSystemError(
    #         "System is not positive",
    #         is_metzler ? nothing : A,
    #         C_positive ? nothing : C
    #     ))
    # end
    
    return true
end

"""
    NonObservableSystemError

Thrown when the system is not observable.

# Fields 
- 'msg::String': Error message
- 'rank_deficit::Int': The system matrix
- 'n::Int': System dimension
"""
struct NonObservableSystemError <: IntervalObserversError
    msg::String
    rank_deficit::Int
    n::Int
end

# function Base.showerror(io::IO, e::NonObservableSystemError)
#     println(io, "NonObservableSystemError: ", e.msg)
# end 


"""
    InvalidInitialBoundsError

Thrown when initial bounds don't satisfy x⁻₀ ≤ x₀ ≤ x⁺₀.
"""
struct InvalidInitialBoundsError <: IntervalObserversError
    violations::Vector{Tuple{Int, Float64,  Float64,  Float64}}
end

""""
    validate_initial_bounds(x0::Vector, xl0::Vector, xu0::Vector)
Throws DimensionMismatchError when initial bounds don't satisfy x⁻₀ ≤ x₀ ≤ x⁺₀.
"""
function validate_initial_bounds(x0::Vector, xl0::Vector, xu0::Vector; kwargs...)
    n = length(x0)

    if length(xl0) != n || length(xu0) != n 
        throw(DimensionMismatchError(
            "initial condition dimension don't match: x0: $n, xl0: $(length(xl0)), xu0: $(length(xu0))"
        ))
    end 
end



struct NotImplemented <: IntervalObserversError
    var::String
end

"""
Customizes the printed message of the exception.
"""
function Base.showerror(io::IO, e::NotImplemented)
    printstyled(io, "NotImplemented"; color=:red, bold=true)
    return print(io, ": ", e.var)
end

"""
    DimensionMismatchError

Thrown when matrix/vector dimensions don't match
"""
struct DimensionMismatchError <: IntervalObserversError
    msg::String
end

function validate_system_dimensions(A::Matrix, C::Vector; kwargs...)
    n, p = size(A)

    if n != p 
        throw(DimensionMismatchError(
            "A must be square matrix, got size $(size(A))"
        ))
    end 

    if length(C) != n 
        throw(DimensionMismatchError(
            "C lenght ($(length(C))) must match A dimensions ($(n), $(n))"
        ))
    end 

    for (key, value) in kwargs
        if key == :state_names
            if length(value) != n
                throw(DimensionMismatchError(
                    "state_names length ($(length(value))) must match A dimensions ($(n), $(n))"
                ))
            end
        end
    end
    return n 
end

"""

"""
function check_observability(M::Matrix, n::Int)
    p = rank(M)

    if p != n
        throw(NonObservableSystemError(
            "The system is not observable: rank(M) = $(p) is different from $n ",
            n - p,
            n
        ))
    end 

    return true
end


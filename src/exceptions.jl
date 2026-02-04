"""
    IntervalObserverError

Abstract type for all interval observer-related errors.
"""

abstract type IntervalObserversError <: Exception end


"""
    NonPositiveSystemError

Thrown when attemting to create an interval observer for non positive system
# Fields 
- 'msg::String': Error message
- 'A::Matrix': The system matrix
- 'C::Vector': The output vector 
"""
struct NonPositiveSystemError <: IntervalObserversError
    msg::String
    A::Union{Nothing, Matrix}
    C::Union{Nothing, Vector}
end

function Base.showerror(io::IO, e::NonPositiveSystemError)
    printstyled(io, "NonPositiveSystemError"; color=:red, bold=true)
    print(io, ": ", e.msg, "\n")
    
    if !isnothing(e.A)
        n = size(e.A, 1)
        non_metzler = Tuple{Int,Int,Float64}[]
        
        # Check for non-Metzler entries
        for j in 1:n, i in 1:n
            if i ≠ j && e.A[i,j] > 0
                push!(non_metzler, (i, j, e.A[i,j]))
            end
        end
        
        if !isempty(non_metzler)
            printstyled(io, "  Non-Metzler entries in A :\n"; color=:yellow)
            for (i, j, val) in non_metzler
                print(io, "    A[$i,$j] = $val > 0\n")
            end
        end
    end
    
    if !isnothing(e.C)
        non_positive_c = findall(x -> x < 0, e.C)
        if !isempty(non_positive_c)
            printstyled(io, "  Non-positive entries in C :\n"; color=:yellow)
            for idx in non_positive_c
                print(io, "    C[$idx] = $(e.C[idx]) < 0\n")
            end
        end
    end
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

function Base.showerror(io::IO, e::NonObservableSystemError)
    println(io, "NonObservableSystemError: ", e.msg)
end 


"""
    InvalidInitialBoundsError

Thrown when initial bounds don't satisfy x⁻₀ ≤ x₀ ≤ x⁺₀.

"""
struct InvalidInitialBoundsError <: IntervalOserversError
    violations::Vector{Tuple{Int, Float64,  Float64,  Float64}}
end

function Base.showerror(io::IO, e::InvalidInitialBoundsError)
    println(io, "InvalidInitialBoundsError: Initial conditions must satisfy x⁻₀ ≤ x₀ ≤ x⁺₀")
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
struct ObservedIntervalProblem{T<:Real}
    t0::T
    tf::T
    n::Int
    state_names::Vector{String}
    A::Matrix{T}
    C::Vector{T}
    xl::Vector{T}
    xu::Vector{T}
    f::Union{Function, Nothing}
    is_linear::Bool
    
    function ObservedIntervalProblem(; t0::Real, 
                                     tf::Real, 
                                     state_names::Vector{String},
                                     A::Matrix{<:Real}, 
                                     C::Vector{<:Real},
                                     xl::Vector{<:Real},
                                     xu::Vector{<:Real},
                                     f::Union{Function, Nothing} = nothing)
        
        # Infer n from state_names length
        n = length(state_names)
        
        # Promote all numeric types to a common type
        T = promote_type(eltype(A), eltype(C), eltype(xl), eltype(xu), typeof(t0), typeof(tf))
        
        # Convert everything to common type
        t0_converted = convert(T, t0)
        tf_converted = convert(T, tf)
        A_converted = convert(Matrix{T}, A)
        C_converted = convert(Vector{T}, C)
        xl_converted = convert(Vector{T}, xl)
        xu_converted = convert(Vector{T}, xu)
        
        # Validation
        @assert tf_converted > t0_converted "Final time must be greater than initial time: tf=$tf_converted, t0=$t0_converted"
        @assert size(A_converted) == (n, n) "A must be $(n)×$(n), got $(size(A_converted))"
        @assert length(C_converted) == n "C must have length $(n), got $(length(C_converted))"
        @assert length(xl_converted) == n "xl must have length $(n), got $(length(xl_converted))"
        @assert length(xu_converted) == n "xu must have length $(n), got $(length(xu_converted))"
        
        # Check interval constraints
        for i in 1:n
            @assert xl_converted[i] < xu_converted[i] "Lower bound must be less than upper bound for state $(state_names[i]): xl[$(i)]=$(xl_converted[i]) >= xu[$(i)]=$(xu_converted[i])"
        end
        
        is_linear = f === nothing
        
        new{T}(t0_converted, tf_converted, n, state_names, A_converted, C_converted, xl_converted, xu_converted, f, is_linear)
    end
end

macro def(expr)
    # Parse the expression to extract time interval syntax
    function parse_def_block(ex)
        if ex isa Expr && ex.head == :block
            new_exprs = []
            for item in ex.args
                if item isa Expr && item.head == :call && item.args[1] == :∈
                    # Skip: t ∈ [t0, tf] - just ignore this line, t0 and tf should already be defined
                    continue
                else
                    push!(new_exprs, item)
                end
            end
            return Expr(:block, new_exprs...)
        end
        return ex
    end
    
    parsed_expr = parse_def_block(expr)
    
    quote
        let
            # Initialize all variables as nothing
            local t0 = nothing
            local tf = nothing
            local z = nothing
            local states = nothing
            local A = nothing
            local C = nothing
            local xl = nothing
            local xu = nothing
            local f = nothing
            
            # Evaluate the parsed block
            $(esc(parsed_expr))
            
            # Handle different syntaxes for time interval
            if t0 === nothing || tf === nothing
                error("Time interval must be specified using 't ∈ [t0, tf]' or 't0 = ..., tf = ...'")
            end
            
            # Handle different syntaxes for states
            state_names = if z !== nothing
                z
            elseif states !== nothing
                states
            else
                error("State variables must be specified using 'z = [...]' or 'states = [...]'")
            end
            
            # Validate required fields
            if A === nothing
                error("Matrix A must be defined")
            end
            if C === nothing
                error("Vector C must be defined")
            end
            if xl === nothing
                error("Lower bounds xl must be defined")
            end
            if xu === nothing
                error("Upper bounds xu must be defined")
            end
            
            # Create and return the OIP
            ObservedIntervalProblem(
                t0 = t0,
                tf = tf,
                state_names = state_names,
                A = A,
                C = C,
                xl = xl,
                xu = xu,
                f = f
            )
        end
    end
end

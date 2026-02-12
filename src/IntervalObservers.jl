"""
IntervalObservers module.

List of exported functions:

"""


module IntervalObservers

using DocStringExtensions

include("exceptions.jl")
include("nonlinear/system.jl")
include("linear/system.jl")
include("linear/simulation.jl")
# include("linear/interval_observer.jl")
include("linear/visualization.jl")
include("IntervalObserversProblem.jl")
include("solve.jl")

export LinearSystem, 
       compute_observability_matrix,
       plot_state_intervals,
       Linear_syst_int_obs_ode!,
       solve, 
       get_state,
       get_upper, 
       get_lower,
       IntervalObserversError,
       validate_system_dimensions, 
       DimensionMismatchError, 
       @def,
       IntervalObserversProblem

# macro def(exp)
#        equations = exp.args

#        oip = Dict(equations => equations)

#        return esc(oip)
# end
end 

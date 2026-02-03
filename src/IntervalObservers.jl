"""
IntervalsObservers module.

List of exported functions:

"""


module IntervalsObservers

using DocStringExtensions

include("nonlinear/system.jl")
include("linear/system.jl")
include("linear/simulation.jl")
# include("linear/interval_observer.jl")
include("linear/visualization.jl")
include("solve.jl")

export LinearSystem, 
       compute_observability_matrix,
       plot_state_intervals,
       Linear_syst_int_obs_ode!,
       solve, 
       get_state,
       get_upper, 
       get_lower
end 

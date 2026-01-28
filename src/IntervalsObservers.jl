"""IntervalsObservers module.

List of exported functions:

"""


module IntervalsObservers

using DocStringExtensions

include("nonlinear/system.jl")
include("linear/system.jl")
include("linear/simulation.jl")
include("linear/interval_observer.jl")
include("linear/visualization.jl")

export LinearSystem, 
       is_observable, 
       compute_observability_matrix,
       observer_gain,
       state_observer,
       simulate_observer,
       observer_info,
       convergence_time,
       IntervalObserver,
       create_interval_observer,
       simulate_interval_observer,
       interval_width,
       check_bounds_validity,
       interval_observer_info,
       plot_states,
       plot_output,
       plot_interval_states,
       plot_interval_width,
       plot_interval_summary


end 

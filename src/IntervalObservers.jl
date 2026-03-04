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
include("nonlinear/simulation.jl")
include("visualization.jl")
include("IntervalObserversProblem.jl")
include("solve.jl")
include("utility.jl")

export LinearSystem, 
       compute_observability_matrix,
       plot_state_intervals,
       plot_nonlinear_state_intervals,
       Linear_syst_int_obs_ode!,
       NonLinear_syst_int_obs_ode!,
       solve, 
       get_state,
       get_upper, 
       get_lower,
       IntervalObserversError,
       NonMetzlerMatrixError,
       NonObservableSystemError,
       InvalidInitialBoundsError,
       DimensionMismatchError,
       NonMonotoneDynamicsError,
       UnobservableMeasurementError,
       InvalidDesiredPolesError,
       validate_system_dimensions, 
       IntervalObserversProblem,
       NonLinearSystem,
       IntervalObserver,
       observable_canonical_form,
       interval_observer_gain,
       evaluate_f,
       get_upper_nonlinear,
       get_lower_nonlinear, 
       positive_interval_gain
end 

# # Nonlinear System Example
#
# This example shows how to use `IntervalObservers.jl` for a nonlinear system
# with interval-bounded uncertainty.
#
# ## Problem Setup
#
# We consider the nonlinear system
#
# ```math
# \dot{x}(t) = Ax(t) + f(t, x(t)),
# ```
#
# with measured output
#
# ```math
# y(t) = Cx(t).
# ```
#
# The nonlinear term is bounded by lower and upper interval functions
# $f_{-}(t,x)$ and $f_{+}(t,x)$.

using IntervalObservers

# ## Parameters
#
# First, define the model parameters.

α₁ = 0.5
m1 = 0.1
α₂ = 0.3
m2 = 0.1
m3 = 0.05
β  = 1.0

# ## Linear Dynamics
#
# Construct the system matrix `A` and output matrix `C`.

β₁ = α₁ + m1
β₂ = α₂ + m2
β₃ = m3

A = [-β₁ 0.0 0.0;
      α₁ -β₂ 0.0;
      0.0  α₂ -β₃]

C = [0.0; 0.0; 1.0]

# The matrix `A` represents the linear dynamics and `C` selects the third state
# as the measured output.

# ## Nonlinear Uncertainty Functions
#
# We define Michaelis-Menten-type nonlinearities.

function f₁(t, a, y)
    a * y / (1.0 + y)
end

function f₃(t, c, y)
    -c * y
end

# Define parameter bounds.

a_max = 0.4
a_min = 0.1
c_max = 0.1
c_min = 0.01

# Define the upper bound function $f_{+}$.

f_plus = [
    (t, y) -> f₁(t, a_max, y),
    (t, y) -> 0.0,
    (t, y) -> f₃(t, c_max, y),
]

# Define the lower bound function $f_{-}$.

f_minus = [
    (t, y) -> 0.5 * f₁(t, a_min, y),
    (t, y) -> 0.0,
    (t, y) -> 0.5 * f₃(t, c_min, y),
]

# ## Create the Nonlinear System

sys = NonLinearSystem(A, C, f_plus, f_minus)
tspan = (0.0, 3.0)

# ## Initial Conditions
#
# Define lower and upper initial interval bounds.

x0_minus = [0.05; 0.1; 0.15]
x0_plus  = [0.15; 0.3; 0.45]

# Define the true initial state.

x0 = [0.1; 0.2; 0.3]

# ## Observer Gain
#
# Compute the observer gain using desired pole locations.

K = positive_interval_gain(sys, desired_poles = [-1.0, -2.0, -3.0])

# ## Solve the Observer Problem

sol_observer = IntervalObservers.solve(
    sys,
    K,
    f_plus,
    f_minus,
    x0_plus,
    x0_minus,
    tspan;
    x0 = x0,
)

# ## Visualization
#
# The following plot shows the true trajectory together with the lower and upper
# interval estimates.

plot_nonlinear_state_intervals(sol_observer, sys)
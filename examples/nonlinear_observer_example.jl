using IntervalObservers

# ## Parameters
# First, define the model parameters.

α₁ = 0.5
α₂ = 0.3

m1 = 0.1
m2 = 0.1
m3 = 0.05


β  = 1.0
b = 20.0

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

# Define parameter bounds.

a_min = 0.1
a_max = 0.6

c_min = 0.01
c_max = 0.2

# ## Nonlinear Uncertainty Functions
#
# We define Michaelis-Menten-type nonlinearities.

function f₁(t, a, y, b)
    a * y / (b + y)
end

function f₃(t, c, y)
    -c * y
end

f_plus = [
    (t, y) -> f₁(t, a_max, y, b),   # upper recruitment
    (t, y) -> 0.0,
    (t, y) -> f₃(t, c_min, y),      # upper harvesting term
]

f_minus = [
    (t, y) -> f₁(t, a_min, y, b),   # lower recruitment
    (t, y) -> 0.0,
    (t, y) -> f₃(t, c_max, y),      # lower harvesting te
]

a_true = 0.35
c_true = 0.08

f_true = [
    (t, y) -> f₁(t, a_true, y, b),
    (t, y) -> 0.0,
    (t, y) -> f₃(t, c_true, y),
]
# ## Create the Nonlinear System

sys = NonLinearSystem(A, C, f_plus, f_minus)
tspan = (0.0, 3.0)

# ## Initial Conditions
#
# Define lower and upper initial interval bounds.


x0_minus = [2.0; 4.0; 6.0]
x0_plus  = [8.0; 16.0; 24.0]
x0       = [5.0; 10.0; 15.0]

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
    f_true = f_true
)
### Plot the Results
plot_nonlinear_state_intervals(sol_observer, sys)
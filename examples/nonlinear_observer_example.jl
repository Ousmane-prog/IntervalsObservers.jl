using IntervalObservers
# using Plots
β₁ = 0.5
β₂ = 0.3
β₃ = 0.2
α₁ = 0.1
α₂ = 0.1

A = [-β₁ 0.0 0.0;
     α₁ -β₂ 0.0;
     0.0 α₂ -β₃]

C = [0.0; 0.0; 1.0]

function f₁(t, a, y)
    return a*y / (1.0 + y)
end

function f₃(t, c, y)
    return c*y / (1.0 + y)
end

# Define f as a Vector of 3 functions (one for each state)
# f = [
#     (t, y) -> f₁(t, α₁, y),
#     (t, y) -> 0.0,
#     (t, y) -> f₃(t, α₂, y)
# ]

# Define f_plus and f_minus as Vectors of functions
# f_plus: upper bound on the nonlinearity (uses true function)
f_plus = [
    (t, y) -> f₁(t, α₁, y),
    (t, y) -> 0.0,
    (t, y) -> f₃(t, α₂, y)
]

# f_minus: lower bound on the nonlinearity (conservative lower bound)
f_minus = [
    (t, y) -> 0.5 * f₁(t, α₁, y),
    (t, y) -> 0.0,
    (t, y) -> 0.5 * f₃(t, α₂, y)
]

# Create the nonlinear system (only f_plus and f_minus are needed)
sys = NonLinearSystem(A, C, f_plus, f_minus)
tspan = (0.0, 3.0)

# Initial state bounds (x0 is not needed)
x0_minus = [0.05; 0.1; 0.15]    
x0_plus  = [0.15; 0.3; 0.45]    

# Use a default gain for now
K = [0.1; 0.1; 0.2]

# Create the interval observer
obs = IntervalObserver(sys, K, f_plus, f_minus)

# Solve the interval observer problem (only pass the bounds, not x0)
sol = IntervalObservers.solve(obs, x0_plus, x0_minus, tspan)

# Extract bounds
xl = get_lower_nonlinear(sol, sys.n)
xu = get_upper_nonlinear(sol, sys.n)

# Plot the nonlinear interval observer results
plt = plot_nonlinear_state_intervals(sol, obs)
display(plt)

println("Nonlinear interval observer solution computed successfully!")
println("State dimension: ", sys.n)
println("Time span: ", tspan)
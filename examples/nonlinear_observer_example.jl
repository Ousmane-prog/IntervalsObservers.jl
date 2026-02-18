using IntervalObservers

α₁ = 0.5
m1 = 0.1
α₂ = 0.3
m2 = 0.1
m3 = 0.05
β  = 1.0

β₁ = α₁ + m1           
β₂ = α₂ + m2          
β₃ = m3           

A = [-β₁ 0.0 0.0;
     α₁ -β₂ 0.0;
     0.0 α₂ -β₃]

C = [0.0; 0.0; 1.0]

function f₁(t, a, y)
    return a*y / (1.0 + y)
end

function f₃(t, c, y)
    return -c*y
end
a_max = 0.4
a_min = 0.1
c_max = 0.1
c_min = 0.01

f_plus = [
    (t, y) -> f₁(t, a_max, y),
    (t, y) -> 0.0,
    (t, y) -> f₃(t, c_max, y)
]

f_minus = [
    (t, y) -> 0.5 * f₁(t, a_min, y),
    (t, y) -> 0.0,
    (t, y) -> 0.5 * f₃(t, c_min, y)
]

sys = NonLinearSystem(A, C, f_plus, f_minus)
tspan = (0.0, 15.0)

x0_minus = [0.05; 0.1; 0.15]    
x0_plus  = [0.15; 0.3; 0.45]    


# K = [0.1; 0.1; 0.2]
K = positive_interval_gain(sys)



obs = IntervalObserver(sys, K, f_plus, f_minus)
sol_observer = IntervalObservers.solve(obs, x0_plus, x0_minus, tspan)

# Plot the observer-only results
plt1 = plot_nonlinear_state_intervals(sol_observer, obs)
display(plt1)

# println("\nObserver-only solution computed successfully!")
# println("State dimension: ", sys.n)
# println("Time span: ", tspan)

# Example 2: Solve WITH true initial state (for comparison with true trajectory)
# x0 = [0.1; 0.2; 0.3]
# sol_with_x0 = IntervalObservers.solve(obs, x0_plus, x0_minus, tspan; x0=x0)

# # Plot with true state
# plt2 = plot_nonlinear_state_intervals(sol_with_x0, obs)
# display(plt2)

# println("\nSolution with true initial state computed successfully!")
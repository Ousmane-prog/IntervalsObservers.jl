using IntervalObservers

# A = [-1.0  0.3;
#       0.0 -0.5]

# C = [1.0, 0.5]

# A = [-2.0  1.0;
#       1.0 -3.0]

# C = [1.0, 0.0]

# x0  = [1.0, 0.4]
# xl0 = [0.6, 0.1]
# xu0 = [1.4, 0.8]

A = [-3.0  1.0  0.5;
      1.0 -4.0  1.0;
      0.5  1.0 -2.0]

C = [1.0, 0.0, 0.0]

x0  = [1.0, 0.5, 0.3]     
xl0 = [0.6, 0.2, 0.0]      
xu0 = [1.4, 0.8, 0.6] 

sys = LinearSystem(A, C)
 
tspan = (0.0, 3.0)

sol = IntervalObservers.solve(
    sys, x0, xl0, xu0,
    tspan
)

# x  = get_state(sol, sys.n)
# xl = get_lower(sol, sys.n)
# xu = get_upper(sol, sys.n)

# (0.0, 10.0)
plot_state_intervals(sol, sys)


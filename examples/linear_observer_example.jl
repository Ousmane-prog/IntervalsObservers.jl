using IntervalObservers

A = [-3.0  1.0  0.5;
      1.0 -4.0  1.0;
      0.5  1.0 -2.0]

C = [1.0, 0.0, 0.0]

x0  = [1.0, 0.5, 0.3]     
xl0 = [0.6, 0.2, 0.0]      
xu0 = [1.4, 0.8, 0.6] 

sys = LinearSystem(A, C)

# desired_poles = [-2.0, -4.0]

# K = positive_interval_gain(sys)
K = positive_interval_gain(sys, desired_poles=[-2.0, -4.0, -6.0])

 
tspan = (0.0, 3.0)

sol = IntervalObservers.solve(
    sys, x0, xl0, xu0,
    K, tspan
)
# redaction 
# (0.0, 10.0)
plot_state_intervals(sol, sys)


# ========================================================================
# ========================================================================

# sys = NonLinearSystem(A, C, f_plus, f_minus)
# tspan = (0.0, 15.0)

# x0_minus = [0.05; 0.1; 0.15]    
# x0_plus  = [0.15; 0.3; 0.45]    


# # K = [0.1; 0.1; 0.2]
# K = positive_interval_gain(sys)



# iop = IntervalObserver(sys, K, f_plus, f_minus)
# sol_observer = IntervalObservers.solve(iop, x0_plus, x0_minus, tspan)

# # Plot the observer-only results
# plt1 = plot_nonlinear_state_intervals(sol_observer, iop)
# display(plt1)






# oip_linear = @def begin
#     t0= 0.0
#     tf = 10.0
#     t ∈ [t0, tf] 
#     z = ["x2", "x2"]
#     A = [0.0 1.0; 
#          -2.0 -0.5]
#     C = [1.0, 0.0]
#     xl = [-2.0, -3.0]
#     xu = [2.0, 3.0]
# end

# (t0, tf) = get_time_interval(oip_linear)
# println("Time: t ∈ [$t0, $tf]")


# @def begin
    
# end
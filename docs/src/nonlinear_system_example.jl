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
    (t, y) -> f₃(t, c_max, y),      # lower harvesting term
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
λ_vals = (-10, -1)
n = size(A, 1)
####################################
######## Utility functions ########
function create_collection(λ_vals::Tuple{Float64, Float64}, n::Int)
    λ_min, λ_max = λ_vals
    @assert λ_min < 0 "λ_min must be negative for stability"
    @assert λ_max < 0 "λ_max must be negative for stability"
    @assert λ_min < λ_max "λ_min must be less than λ_max"

    return collect(λ_min, λ_max, length=n)
end

function generate_poles_geometric(λ::Float64, n::Int; δ=0.5)
    @assert λ < 0 "λ must be negative for stability"
    @assert 0 < δ < 1 "δ must be in (0, 1) for distinct poles"

    return [λ_k_plus_1 = λₖ*(1-δ) for k in 0:(n-1)]
    
end

function generate_poles(λ_vals::Tuple{Float64, Float64}, n::Int)
    pole_collection = create_collection(λ_vals, n)
    poles_list = []
    for λ_vect in pole_collection
        poles = generate_poles_geometric(λ_vect, n)
        push!(poles_list, poles)
    end
    return poles_list
end
######################################
######################################

desired_poles = generate_poles(λ_vals, n)

for poles in desired_poles
    println("Desired poles: ", poles)
    K = positive_interval_gain(sys, desired_poles = poles)
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
    plot_nonlinear_state_intervals(sol_observer, sys)
end
K = positive_interval_gain(sys, desired_poles = [-1.0, -2.0, -3.0])

# ## Solve the Observer Problem

# sol_observer = IntervalObservers.solve(
#     sys,
#     K,
#     f_plus,
#     f_minus,
#     x0_plus,
#     x0_minus,
#     tspan;
#     x0 = x0,
#     f_true = f_true
# )
# ### Plot the Results
# plot_nonlinear_state_intervals(sol_observer, sys)
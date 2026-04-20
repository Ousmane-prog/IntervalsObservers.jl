using DifferentialEquations
using Plots

# Parameters
α1 = 0.5
α2 = 0.3
m1 = 0.1
m2 = 0.1
m3 = 0.05
b  = 20.0

# Time-varying or constant biological inputs
a_true_val(t) = 0.35              # recruitment intensity
c_true_val(t) = 0.08              # harvesting effort

# Beverton-Holt recruitment
recruitment(t, x3) = a_true_val(t) * x3 / (b + x3)

# Fish population model
function fish_population!(dx, x, p, t)
    x1, x2, x3 = x

    dx[1] = -(α1 + m1) * x1 + recruitment(t, x3)
    dx[2] =  α1 * x1 - (α2 + m2) * x2
    dx[3] =  α2 * x2 - m3 * x3 - c_true_val(t) * x3
end

# Initial conditionjulia

x0 = [5.0, 10.0, 15.0]
tspan = (0.0, 20.0)

prob = ODEProblem(fish_population!, x0, tspan)
sol = DifferentialEquations.solve(prob, Tsit5())

plot(sol, xlabel="time", ylabel="population", label=["x₁ larvae" "x₂ juveniles" "x₃ adults"])
# savefig("sol.png")
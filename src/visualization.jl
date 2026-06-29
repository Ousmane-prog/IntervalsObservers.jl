using Plots
import Plots: plot, plot!
using RecipesBase

"""
    IntervalObserverSolution

Stores and represents the solution from an interval observer simulation.

Captures the time trajectory and state estimates (upper/lower bounds, and optionally
the true state) from an interval observer run.

# Fields
- `t`: Time vector
- `u`: State solution as a vector of state vectors
- `label::String`: Description of the observer (e.g., "linear observer")
- `show_true::Bool`: Whether to display the true state trajectory in plots

# Constructor
```julia
IntervalObserverSolution(t, u; label="", show_true=true)
```

# Plotting
The solution can be plotted using Julia's `plot()` function:
```julia
sol = solve(...)  # Returns IntervalObserverSolution
plot(sol)
```

Automatically generates subplots for each state showing:
- Shaded interval between lower and upper bounds
- Dashed lines for upper and lower bounds
- True state trajectory (if available)
"""
struct IntervalObserverSolution
    t
    u
    label::String
    show_true::Bool
end

IntervalObserverSolution(t, u; label="", show_true=true) =
    IntervalObserverSolution(t, u, String(label), show_true)

"""
    plot(sol::IntervalObserverSolution, sys; kwargs...)

Plot the interval observer solution.

Creates a figure with subplots for each state dimension showing:
- Shaded region between lower and upper bounds
- Dashed lines indicating bound trajectories
- True state (if available)

# Arguments
- `sol::IntervalObserverSolution`: The solution to plot
- `sys`: System (optional, for display purposes)
- `kwargs`: Additional plotting options passed to Plots.jl

# Returns
- `Plot`: Plots.jl figure object
"""
function plot(sol::IntervalObserverSolution, sys; kwargs...)
    return plot(sol; kwargs...)
end

"""
    plot!(plt::Plots.Plot, sol::IntervalObserverSolution, sys; kwargs...)

Add an interval observer solution to an existing plot.

# Arguments
- `plt::Plots.Plot`: Existing plot object
- `sol::IntervalObserverSolution`: Solution to add
- `sys`: System (optional)
- `kwargs`: Additional plotting options

# Returns
- `Plot`: Modified plot object
"""
function plot!(plt::Plots.Plot, sol::IntervalObserverSolution, sys; kwargs...)
    return plot!(plt, sol; kwargs...)
end

@recipe function f(sol::IntervalObserverSolution)
    t = sol.t
    Z = reduce(hcat, sol.u)

    num_states = size(Z, 1)
    has_true = iszero(num_states % 3)

    if has_true
        n = num_states ÷ 3
        x  = Z[1:n, :]
        xu = Z[n+1:2n, :]
        xl = Z[2n+1:3n, :]
    else
        n = num_states ÷ 2
        x  = nothing
        xu = Z[1:n, :]
        xl = Z[n+1:2n, :]
    end

    layout := (n, 1)
    size := (800, 250n)
    xlabel := "time (s)"
    legend := false

    for i in 1:n
        lower = min.(xl[i, :], xu[i, :])
        upper = max.(xl[i, :], xu[i, :])

        # shaded interval
        @series begin
            subplot := i
            fillrange := upper
            fillalpha := 0.12
            label := ""
            t, lower
        end

        # upper
        @series begin
            subplot := i
            linestyle := :dash
            linewidth := 1.5
            # label := "upper ($(sol.label))"
            label := nothing
            t, upper
        end

        # lower
        @series begin
            subplot := i
            linestyle := :dash
            linewidth := 1.5
            # label := "lower ($(sol.label))"
            label := nothing
            t, lower
        end

        # true (only once)
        if has_true && sol.show_true
            @series begin
                subplot := i
                linewidth := 2
                color := :black
                label := "true"
                t, x[i, :]
            end
        end

        ylabel := "state $i"
    end
end
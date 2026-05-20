using Plots
import Plots: plot, plot!
using RecipesBase

struct IntervalObserverSolution
    t
    u
    label::String
    show_true::Bool
end

IntervalObserverSolution(t, u; label="", show_true=true) =
    IntervalObserverSolution(t, u, String(label), show_true)

function plot(sol::IntervalObserverSolution, sys; kwargs...)
    return plot(sol; kwargs...)
end

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
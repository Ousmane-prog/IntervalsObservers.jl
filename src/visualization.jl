using Plots
import Plots: plot, plot!

# # Initialize GR backend with appropriate settings for various environments
# try
#     gr()
#     Plots.GR.inline(true)
#     # Set to use agg backend for safer non-interactive rendering
#     Plots.GR.batchmode(false)
# catch
#     # If GR initialization fails, continue with whatever backend is available
# end

function plot_state_intervals(sol, sys)
    n = sys.n
    t = sol.t

    x = get_state(sol, n)
    xl = get_lower(sol, n)
    xu = get_upper(sol, n)
    println("Plotting state intervals for linear observer.")
    plt = plot(layout = (n, 1), size=(800, 250*n))

    for i in 1:n
        plot!(
            plt[i], 
            t, xu[i, :],
            label = "x⁺_$i",
            ls = :dash,
            color = :red
        )
        plot!(
            plt[i], 
            t, x[i, :],
            label = "x_$i",
            lw = 2,
            color = :black
        )

        plot!(
            plt[i], 
            t, xl[i, :],
            label = "x⁻_$i",
            ls = :dash,
            color = :blue
        )

        plot!(
            plt[i], 
            t, xu[i, :],
            fillrange = xl[i,:],
            fillalpha = 0.15,
            color = :lightblue,
            label = nothing
        )
        ylabel!(plt[i], "x_$i")
    end 
    xlabel!(plt[end], "time")
    return plt
end


"""
    plot_nonlinear_state_intervals(sol, sys)

Plot state interval bounds for a nonlinear interval observer solution.

# Arguments
- `sol`: Solution from solver
- `sys::Union{LinearSystem, NonLinearSystem}`: System information

# Returns
- Plot object with subplots for each state dimension
"""
# function plot_nonlinear_state_intervals(sol, sys)

#     n = sys.n
#     t = sol.t

#     println("Plotting state intervals for nonlinear observer.")

#     num_states = size(sol,1)
#     track_true_state = (num_states == 3n)

#     # println("Number of states in solution: ", num_states)
#     # println("Tracking true state: ", track_true_state)

#     if track_true_state
#         x  = get_state(sol,n)
#         xl = get_lower(sol,n)
#         xu = get_upper(sol,n)
#     else
#         xl = get_lower_nonlinear(sol,n)
#         xu = get_upper_nonlinear(sol,n)
#     end

#     plt = plot(layout=(n,1), size=(800,250n))

#     for i in 1:n

#         lower = min.(xl[i,:], xu[i,:])
#         upper = max.(xl[i,:], xu[i,:])

#         # shaded interval
#         plot!(
#             plt[i],
#             t, lower,
#             fillrange = upper,
#             fillalpha = 0.2,
#             color = :lightblue,
#             label = nothing
#         )

#         # upper bound
#         plot!(
#             plt[i],
#             t, upper,
#             ls = :dash,
#             color = :red,
#             lw = 1.5,
#             label = "upper"
#         )

#         # lower bound
#         plot!(
#             plt[i],
#             t, lower,
#             ls = :dash,
#             color = :blue,
#             lw = 1.5,
#             label = "lower"
#         )

#         if track_true_state
#             plot!(
#                 plt[i],
#                 t, x[i,:],
#                 lw = 2,
#                 color = :black,
#                 label = "true"
#             )
#         end

#         ylabel!(plt[i], "state $i")
#     end

#     xlabel!(plt[end], "time (s)")

#     title_text = track_true_state ?
#         "Nonlinear Interval Observer" :
#         "Nonlinear Interval Observer Bounds"
    
#     # title_text *= "\n[Plot shown in transformed coordinates]"
    
#     title!(plt, title_text, subplot=1)

#     return plt
# end

function plot_nonlinear_state_intervals(sol, sys)
    n = sys.n
    t = sol.t
    Z = reduce(hcat, sol.u)

    println("Plotting state intervals for nonlinear observer.")

    num_states = size(Z, 1)
    track_true_state = (num_states == 3n)

    if track_true_state
        x  = get_state(Z, n)
        xl = get_lower(Z, n)
        xu = get_upper(Z, n)
    else
        xl = get_lower_nonlinear(Z, n)
        xu = get_upper_nonlinear(Z, n)
        x = nothing
    end

    plt = plot(layout=(n,1), size=(800,250n), legend=false)

    for i in 1:n
        lower = min.(xl[i,:], xu[i,:])
        upper = max.(xl[i,:], xu[i,:])

        plot!(
            plt[i],
            t, lower,
            fillrange=upper,
            fillalpha=0.2,
            color=:lightblue,
            label=nothing
        )

        plot!(plt[i], t, upper, ls=:dash, color=:red,  lw=1.5, label="upper")
        plot!(plt[i], t, lower, ls=:dash, color=:blue, lw=1.5, label="lower")

        if track_true_state
            plot!(plt[i], t, x[i,:], lw=2, color=:black, label="true")
        end

        ylabel!(plt[i], "state $i")
    end

    xlabel!(plt[end], "time (s)")
    title!(plt, track_true_state ? "Nonlinear Interval Observer" :
                                   "Nonlinear Interval Observer Bounds",
           subplot=1)

    return plt
end


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
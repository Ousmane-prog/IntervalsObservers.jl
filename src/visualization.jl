using Plots

function plot_state_intervals(sol, sys)
    n = sys.n
    t = sol.t

    x = get_state(sol, n)
    xl = get_lower(sol, n)
    xu = get_upper(sol, n)

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
    plot_nonlinear_state_intervals(sol, obs::IntervalObserver)

Plot state interval bounds for a nonlinear interval observer solution.

# Arguments
- `sol`: Solution from solver
- `obs::IntervalObserver`: Interval observer containing the system information

# Returns
- Plot object with subplots for each state dimension
"""
function plot_nonlinear_state_intervals(sol, obs::IntervalObserver)
    sys = obs.sys
    n = sys.n
    t = sol.t

    xl = get_lower_nonlinear(sol, n)
    xu = get_upper_nonlinear(sol, n)

    plt = plot(layout = (n, 1), size=(800, 250*n))

    for i in 1:n
        plot!(
            plt[i], 
            t, xu[i, :],
            label = "x⁺_$i",
            ls = :dash,
            color = :red,
            lw = 1.5
        )

        plot!(
            plt[i], 
            t, xl[i, :],
            label = "x⁻_$i",
            ls = :dash,
            color = :blue,
            lw = 1.5
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
        # grid!(plt[i], true, alpha=0.3)
    end 
    xlabel!(plt[end], "time (s)")
    title!(plt, "Nonlinear Interval Observer State Bounds", subplot=1)
    return plt
end
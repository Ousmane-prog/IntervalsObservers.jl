using Plots

# Initialize GR backend with appropriate settings for various environments
try
    gr()
    Plots.GR.inline(true)
    # Set to use agg backend for safer non-interactive rendering
    Plots.GR.batchmode(false)
catch
    # If GR initialization fails, continue with whatever backend is available
end

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
function plot_nonlinear_state_intervals(sol, sys)

    n = sys.n
    t = sol.t

    println("Plotting state intervals for nonlinear observer.")

    num_states = size(sol,1)
    track_true_state = (num_states == 3n)

    # println("Number of states in solution: ", num_states)
    # println("Tracking true state: ", track_true_state)

    if track_true_state
        x  = get_state(sol,n)
        xl = get_lower(sol,n)
        xu = get_upper(sol,n)
    else
        xl = get_lower_nonlinear(sol,n)
        xu = get_upper_nonlinear(sol,n)
    end

    plt = plot(layout=(n,1), size=(800,250n))

    for i in 1:n

        lower = min.(xl[i,:], xu[i,:])
        upper = max.(xl[i,:], xu[i,:])

        # shaded interval
        plot!(
            plt[i],
            t, lower,
            fillrange = upper,
            fillalpha = 0.2,
            color = :lightblue,
            label = nothing
        )

        # upper bound
        plot!(
            plt[i],
            t, upper,
            ls = :dash,
            color = :red,
            lw = 1.5,
            label = "upper"
        )

        # lower bound
        plot!(
            plt[i],
            t, lower,
            ls = :dash,
            color = :blue,
            lw = 1.5,
            label = "lower"
        )

        if track_true_state
            plot!(
                plt[i],
                t, x[i,:],
                lw = 2,
                color = :black,
                label = "true"
            )
        end

        ylabel!(plt[i], "state $i")
    end

    xlabel!(plt[end], "time (s)")

    title_text = track_true_state ?
        "Nonlinear Interval Observer" :
        "Nonlinear Interval Observer Bounds"
    
    # title_text *= "\n[Plot shown in transformed coordinates]"
    
    title!(plt, title_text, subplot=1)

    return plt
end
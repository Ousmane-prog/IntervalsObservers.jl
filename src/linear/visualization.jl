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
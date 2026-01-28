using Plots

"""Helper function to create subscript numbers"""
function subscript(i::Int)
    subscripts = Dict('0'=>'₀', '1'=>'₁', '2'=>'₂', '3'=>'₃', '4'=>'₄',
                     '5'=>'₅', '6'=>'₆', '7'=>'₇', '8'=>'₈', '9'=>'₉')
    return join([get(subscripts, c, c) for c in string(i)])
end

# ============================================================
# Standard Observer Visualization
# ============================================================

function plot_states(result; state_names=nothing, kwargs...)
    n = size(result.x_true, 1)
    
    if state_names === nothing
        state_names = ["x" * subscript(i) for i in 1:n]
    end
    
    plots = []
    for i in 1:n
        p = plot(result.t, result.x_true[i, :], 
                label="True $(state_names[i])",
                linewidth=2,
                linestyle=:solid,
                xlabel="Time [s]",
                ylabel=state_names[i];
                kwargs...)
        plot!(p, result.t, result.x_obs[i, :],
             label="Est. $(state_names[i])",
             linewidth=2,
             linestyle=:dash)
        push!(plots, p)
    end
    
    return plot(plots..., layout=(n, 1), size=(800, 300*n))
end

function plot_output(result; kwargs...)
    p = plot(result.t, result.y,
            label="y(t) = C*x",
            linewidth=2,
            xlabel="Time [s]",
            ylabel="Output";
            kwargs...)
    return p
end

# ============================================================
# Interval Observer Visualization
# ============================================================

"""
    plot_interval_states(result; state_names=nothing, kwargs...)

Plot true states with interval bounds (upper and lower observers).
"""
function plot_interval_states(result; state_names=nothing, kwargs...)
    n = size(result.x_true, 1)
    
    if state_names === nothing
        state_names = ["x" * subscript(i) for i in 1:n]
    end
    
    plots = []
    for i in 1:n
        p = plot(result.t, result.x_upper[i, :], 
                fillrange=result.x_lower[i, :],
                fillalpha=0.3,
                fillcolor=:lightblue,
                label="Interval bounds",
                linecolor=:blue,
                linestyle=:dash,
                xlabel="Time [s]",
                ylabel=state_names[i];
                kwargs...)
        
        plot!(p, result.t, result.x_lower[i, :],
             linecolor=:blue,
             linestyle=:dash,
             label="")
        
        plot!(p, result.t, result.x_true[i, :],
             label="True $(state_names[i])",
             linewidth=2,
             linecolor=:red)
        
        push!(plots, p)
    end
    
    return plot(plots..., layout=(n, 1), size=(800, 300*n))
end

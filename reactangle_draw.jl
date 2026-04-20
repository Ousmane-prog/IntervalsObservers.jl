using Plots

function generate_extreme_points(x_minus::AbstractVector, x_plus::AbstractVector)
    @assert length(x_minus) == length(x_plus)

    function gen(i)
        if i > 1
            E1 = gen(i - 1)
            nb_col1 = size(E1, 2)

            top = vcat(E1, x_minus[i] .* ones(1, nb_col1))
            bottom = vcat(E1, x_plus[i] .* ones(1, nb_col1))

            hcat(top, bottom)
        else
            reshape([x_minus[1], x_plus[1]], 1, 2)
        end
    end

    gen(length(x_minus))
end

function compute_bounds_in_new_basis(
    M::AbstractMatrix,
    x_minus::Union{AbstractVector, AbstractMatrix},
    x_plus::Union{AbstractVector, AbstractMatrix},
    x::Union{AbstractVector, AbstractMatrix, Nothing} = nothing
    )
    # Single box
    if x_minus isa AbstractVector && x_plus isa AbstractVector
        E = generate_extreme_points(x_minus, x_plus)
        EZ = M * E

        z_minus = vec(minimum(EZ, dims=2))
        z_plus  = vec(maximum(EZ, dims=2))
        z       = isnothing(x) ? nothing : M * x

        return z_minus, z_plus, z
    end

    # Several boxes stored columnwise
    @assert x_minus isa AbstractMatrix && x_plus isa AbstractMatrix
    @assert size(x_minus) == size(x_plus)
    @assert isnothing(x) || size(x) == size(x_minus)

    n_out = size(M, 1)
    n_cols = size(x_minus, 2)

    z_minus = Matrix{eltype(x_minus)}(undef, n_out, n_cols)
    z_plus  = Matrix{eltype(x_plus)}(undef, n_out, n_cols)
    z       = isnothing(x) ? nothing : Matrix{eltype(x)}(undef, n_out, n_cols)

    for k in 1:n_cols
        E = generate_extreme_points(x_minus[:, k], x_plus[:, k])
        EZ = M * E

        z_minus[:, k] = vec(minimum(EZ, dims=2))
        z_plus[:, k]  = vec(maximum(EZ, dims=2))

        if !isnothing(x)
            z[:, k] = M * x[:, k]
        end
    end

    return z_minus, z_plus, z
end

function draw_polytope!(
    plt,
    points::AbstractMatrix,
    reference_vertices::AbstractMatrix;
    label_str = "",
    linestyle = :solid,
    markersize = 4, 
    color = :blue
)
    @assert size(points, 1) == 3
    @assert size(reference_vertices, 1) == 3
    @assert size(points, 2) == size(reference_vertices, 2)

    scatter!(
        plt,
        points[1, :], points[2, :], points[3, :],
        label = label_str,
        markersize = markersize,
    )

    m = size(points, 2)
    for i in 1:m
        for j in (i + 1):m
            if sum(reference_vertices[:, i] .!= reference_vertices[:, j]) == 1
                plot!(
                    plt,
                    [points[1, i], points[1, j]],
                    [points[2, i], points[2, j]],
                    [points[3, i], points[3, j]],
                    label = false,
                    linestyle = linestyle,
                    color = color
                )
            end
        end
    end

    plt
end

function plot_box_pipeline(
    M::AbstractMatrix,
    x_minus::AbstractVector,
    x_plus::AbstractVector
)
    # @assert length(x_minus) == 3 "This plotting function is for 3D boxes."
    # @assert length(x_plus) == 3 "This plotting function is for 3D boxes."
    # @assert size(M) == (3, 3) "This plotting function expects a 3×3 matrix."

    # Original box in x-space
    E = generate_extreme_points(x_minus, x_plus)

    # Exact image in z-space
    EZ = M * E

    # Bounding box in z-space
    z_minus, z_plus, _ = compute_bounds_in_new_basis(M, x_minus, x_plus)
    @info "z_minus: $z_minus"
    @info "z_plus: $z_plus"
    EZ_box = generate_extreme_points(z_minus, z_plus)

    # Exact inverse image of the z-box vertices
    M_inv = inv(M)
    # EX_back_exact = M_inv * EZ_box
    EX_back_exact = M_inv * EZ_box

    # Bounding box in x-space of that inverse image
    # x_back_minus = vec(minimum(EX_back_exact, dims = 2))
    # x_back_plus  = vec(maximum(EX_back_exact, dims = 2))
    x_back_minus, x_back_plus, _ = compute_bounds_in_new_basis(M_inv, z_minus, z_plus)
    @info "x_back_minus: $x_back_minus"
    @info "x_back_plus: $x_back_plus"
    EX_back_box = generate_extreme_points(x_back_minus, x_back_plus)

    # x-space plot
    plt_x = plot(
        xlabel = "x1",
        ylabel = "x2",
        zlabel = "x3",
        title = "x-space",
        legend = :topright
    )

    draw_polytope!(
        plt_x, E, E;
        label_str = "original box",
        linestyle = :solid,
        markersize = 5,
        color = :blue
    )


    plt_x_inv = plot(
        xlabel = "x1",
        ylabel = "x2",
        zlabel = "x3",
        title = "x-space",
        legend = :topright
    )

    draw_polytope!(
        plt_x_inv, EX_back_exact, EZ_box;
        label_str = "inverse image exact",
        linestyle = :dash,
        markersize = 4,
        color = :red
    )
    

    plt_x_inv_box = plot(
        xlabel = "x1",
        ylabel = "x2",
        zlabel = "x3",
        title = "x-space",
        legend = :topright
    )

    draw_polytope!(
        plt_x_inv_box, EX_back_exact, EZ_box;
        label_str = "inverse image exact",
        linestyle = :dash,
        markersize = 4,
        color = :red
    )

    draw_polytope!(
        plt_x_inv_box, EX_back_box, EX_back_box;
        label_str = "inverse image box",
        linestyle = :dot,
        markersize = 4,
        color = :green

    )

    # z-space plot
    plt_z = plot(
        xlabel = "z1",
        ylabel = "z2",
        zlabel = "z3",
        title = "z-space",
        legend = :topright
    )

    draw_polytope!(
        plt_z, EZ, E;
        label_str = "exact image",
        linestyle = :solid,
        markersize = 5,
        color = :orange
    )

    draw_polytope!(
        plt_z, EZ_box, EZ_box;
        label_str = "image box",
        linestyle = :dash,
        markersize = 4,
        color = :purple
    )

    plot(plt_x, plt_z, plt_x_inv, plt_x_inv_box, layout = (2, 2), size = (1800, 1200))
end

# Example
# x_minus = [0.0, 1.0, 2.0]
# x_plus  = [1.0, 2.0, 3.0]
x_plus  = [8.0; 16.0; 24.0]
# x0       = [5.0; 10.0; 15.0]
x_minus = [2.0; 4.0; 6.0]

# M = [1.0 2.0 0.0;
#      0.0 1.0 1.0;
#      1.0 0.0 1.0]
M = [0.5779489596841483 -2.774155006483916 24.04267672286065; 
    -1.7872045210327752 5.004172658891751 -26.68892084742255; 
    1.956533925082833 -1.565227140066275 3.1304542801325623]

plot_box_pipeline(M, x_minus, x_plus)

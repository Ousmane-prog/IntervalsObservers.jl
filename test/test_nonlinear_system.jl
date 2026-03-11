@testset "NonLinearSystem Construction and Validation" begin

    @testset "NonLinearSystem Creation with Metzler Check" begin
        A = [-2.0 1.0;
             0.5 -3.0]
        C = [1.0, 0.0]
        
        f_plus = [(t, y) -> 0.1*y, (t, y) -> 0.05*y]
        f_minus = [(t, y) -> 0.05*y, (t, y) -> 0.02*y]
        
        sys = NonLinearSystem(A, C, f_plus, f_minus)
        
        @test sys.n == 2
        @test sys.A == A
        @test sys.C == C
        @test length(sys.f_plus) == 2
        @test length(sys.f_minus) == 2
    end

    @testset "NonLinearSystem with Non-Metzler Matrix" begin
        A = [-2.0 -1.0;  # Non-Metzler (has negative off-diagonal)
             0.5 -3.0]
        C = [1.0, 0.0]
        
        f_plus = [(t, y) -> 0.1*y, (t, y) -> 0.05*y]
        f_minus = [(t, y) -> 0.05*y, (t, y) -> 0.02*y]
        
        sys = NonLinearSystem(A, C, f_plus, f_minus; check_metzler=false)
        
        @test sys.n == 2
        @test sys.A == A
    end

    @testset "Dimension Mismatch in NonLinearSystem" begin
        A = [-2.0 1.0;
             0.5 -3.0]
        C = [1.0, 0.0, 0.5]  # Wrong dimension
        
        f_plus = [(t, y) -> 0.1*y, (t, y) -> 0.05*y]
        f_minus = [(t, y) -> 0.05*y, (t, y) -> 0.02*y]
        
        @test_throws DimensionMismatchError NonLinearSystem(A, C, f_plus, f_minus)
    end

    @testset "Nonlinear Functions Mismatch" begin
        A = [-2.0 1.0;
             0.5 -3.0]
        C = [1.0, 0.0]
        
        f_plus = [(t, y) -> 0.1*y, (t, y) -> 0.05*y]  # 2 functions
        f_minus = [(t, y) -> 0.05*y]  # Only 1 function
        
        @test_throws DimensionMismatchError NonLinearSystem(A, C, f_plus, f_minus)
    end
end


@testset "Nonlinear Interval Observer Solver" begin

    @testset "Basic Nonlinear Observer Monotone System" begin
        α₁ = 0.5
        m1 = 0.1
        α₂ = 0.3
        m2 = 0.1
        m3 = 0.05

        β₁ = α₁ + m1
        β₂ = α₂ + m2
        β₃ = m3

        A = [-β₁  0.0  0.0;
              α₁ -β₂  0.0;
              0.0  α₂ -β₃]

        C = [0.0; 0.0; 1.0]

        function f₁(t, a, y)
            return a*y / (1.0 + y)
        end

        function f₃(t, c, y)
            return -c*y
        end

        a_max = 0.4
        a_min = 0.1
        c_max = 0.1
        c_min = 0.01

        f_plus = [
            (t, y) -> f₁(t, a_max, y),
            (t, y) -> 0.0,
            (t, y) -> f₃(t, c_max, y)
        ]

        f_minus = [
            (t, y) -> 0.5 * f₁(t, a_min, y),
            (t, y) -> 0.0,
            (t, y) -> 0.5 * f₃(t, c_min, y)
        ]

        sys = NonLinearSystem(A, C, f_plus, f_minus)
        tspan = (0.0, 5.0)

        x0_minus = [0.05; 0.1; 0.15]
        x0_plus = [0.15; 0.3; 0.45]
        x0 = [0.1; 0.2; 0.3]

        K = positive_interval_gain(sys, desired_poles = [-1.0, -2.0, -3.0])

        sol = IntervalObservers.solve(sys, K, f_plus, f_minus, x0_plus, x0_minus, tspan; x0=x0)

        @test sol.t[1] == tspan[1]
        @test sol.t[end] ≈ tspan[2]
        @test size(sol.u, 1) == 3  # x, x_plus, x_minus tracked
    end

    @testset "Nonlinear Observer Without True State Tracking" begin
        A = [-2.0 1.0 0.0;
             0.5 -3.0 1.0;
             0.0 0.5 -2.0]

        C = [1.0, 0.0, 0.0]

        f_plus = [(t, y) -> 0.1*y, (t, y) -> 0.0, (t, y) -> 0.05*y]
        f_minus = [(t, y) -> 0.05*y, (t, y) -> 0.0, (t, y) -> 0.02*y]

        sys = NonLinearSystem(A, C, f_plus, f_minus; check_metzler=false)
        tspan = (0.0, 2.0)

        x0_minus = [0.1, 0.2, 0.1]
        x0_plus = [0.3, 0.5, 0.3]

        K = [0.1; 0.1; 0.2]

        sol = IntervalObservers.solve(sys, K, f_plus, f_minus, x0_plus, x0_minus, tspan)

        @test sol.t[1] == tspan[1]
        @test sol.t[end] ≈ tspan[2]
        @test size(sol.u, 1) == 6  # x_plus and x_minus only (no true state)
    end

#     @testset "Initial Bounds Validation" begin
#         A = [-2.0 1.0;
#              0.5 -3.0]
#         C = [1.0, 0.0]

#         f_plus = [(t, y) -> 0.1*y, (t, y) -> 0.05*y]
#         f_minus = [(t, y) -> 0.05*y, (t, y) -> 0.02*y]

#         sys = NonLinearSystem(A, C, f_plus, f_minus; check_metzler=false)
#         tspan = (0.0, 1.0)

#         x0_minus = [0.3, 0.5]  # Lower bound > upper bound
#         x0_plus = [0.1, 0.2]

#         @test_throws InvalidInitialBoundsError IntervalObservers.solve(
#             sys, [0.1; 0.1], f_plus, f_minus, x0_plus, x0_minus, tspan
#         )
#     end
# end


# @testset "Nonlinear Extraction Functions" begin

#     @testset "get_upper_nonlinear and get_lower_nonlinear" begin
#         A = [-2.0 1.0;
#              0.5 -3.0]
#         C = [1.0, 0.0]

#         f_plus = [(t, y) -> 0.1*y, (t, y) -> 0.05*y]
#         f_minus = [(t, y) -> 0.05*y, (t, y) -> 0.02*y]

#         sys = NonLinearSystem(A, C, f_plus, f_minus; check_metzler=false)
#         tspan = (0.0, 1.0)

#         x0_minus = [0.1, 0.2]
#         x0_plus = [0.3, 0.5]

#         K = [0.1; 0.1]

#         sol = IntervalObservers.solve(sys, K, f_plus, f_minus, x0_plus, x0_minus, tspan)

#         n = sys.n
#         xu = get_upper_nonlinear(sol, n)
#         xl = get_lower_nonlinear(sol, n)

#         @test size(xu, 1) == n
#         @test size(xl, 1) == n
#         @test size(xu, 2) == length(sol.t)
#         @test size(xl, 2) == length(sol.t)

#         # Check initial conditions
#         @test xu[:, 1] == x0_plus
#         @test xl[:, 1] == x0_minus
#     end
# end


# @testset "Nonlinear Observer Gain Computation" begin

#     @testset "positive_interval_gain for Nonlinear Systems" begin
#         α₁ = 0.5
#         m1 = 0.1
#         α₂ = 0.3
#         m2 = 0.1
#         m3 = 0.05

#         β₁ = α₁ + m1
#         β₂ = α₂ + m2
#         β₃ = m3

#         A = [-β₁  0.0  0.0;
#               α₁ -β₂  0.0;
#               0.0  α₂ -β₃]

#         C = [0.0; 0.0; 1.0]

#         f_plus = [(t, y) -> 0.1*y, (t, y) -> 0.0, (t, y) -> 0.05*y]
#         f_minus = [(t, y) -> 0.05*y, (t, y) -> 0.0, (t, y) -> 0.02*y]

#         sys = NonLinearSystem(A, C, f_plus, f_minus)

#         # without desired poles
#         K1 = positive_interval_gain(sys)
#         @test length(K1) == 3
#         @test all(K1 .>= 0)  # All gains should be positive

#         # with desired poles
#         K2 = positive_interval_gain(sys, desired_poles = [-1.0, -2.0, -3.0])
#         @test length(K2) == 3
#         @test all(K2 .>= 0)
#     end

#     @testset "Gain for Non-Observable System" begin
#         A = [-2.0 0.0;
#              0.0 -3.0]
#         C = [0.0, 1.0]  # Cannot observe first state

#         f_plus = [(t, y) -> 0.1*y, (t, y) -> 0.0]
#         f_minus = [(t, y) -> 0.05*y, (t, y) -> 0.0]

#         @test_throws NonObservableSystemError NonLinearSystem(A, C, f_plus, f_minus)
#     end
# end


# @testset "Nonlinear System Visualization" begin

#     @testset "plot_nonlinear_state_intervals" begin
#         A = [-2.0 1.0 0.0;
#              0.5 -3.0 1.0;
#              0.0 0.5 -2.0]

#         C = [1.0, 0.0, 0.0]

#         f_plus = [(t, y) -> 0.1*y, (t, y) -> 0.0, (t, y) -> 0.05*y]
#         f_minus = [(t, y) -> 0.05*y, (t, y) -> 0.0, (t, y) -> 0.02*y]

#         sys = NonLinearSystem(A, C, f_plus, f_minus; check_metzler=false)
#         tspan = (0.0, 1.0)

#         x0_minus = [0.1, 0.2, 0.1]
#         x0_plus = [0.3, 0.5, 0.3]
#         x0 = [0.2, 0.35, 0.2]

#         K = [0.1; 0.15; 0.1]

#         sol = IntervalObservers.solve(sys, K, f_plus, f_minus, x0_plus, x0_minus, tspan; x0=x0)

#         # Test that plotting doesn't error and returns a plot object
#         plt = plot_nonlinear_state_intervals(sol, sys)
        
#         @test plt !== nothing
#         @test haskey(plt.subplots, 1)  # At least one subplot exists
#     end
# end


# @testset "Nonlinear Observer Properties" begin

#     @testset "Containment Property" begin
#         A = [-2.0 1.0;
#              0.5 -3.0]
#         C = [1.0, 0.0]

#         f_plus = [(t, y) -> 0.1*y, (t, y) -> 0.05*y]
#         f_minus = [(t, y) -> 0.05*y, (t, y) -> 0.02*y]

#         sys = NonLinearSystem(A, C, f_plus, f_minus; check_metzler=false)
#         tspan = (0.0, 2.0)

#         x0_minus = [0.1, 0.2]
#         x0_plus = [0.3, 0.5]
#         x0 = [0.2, 0.35]

#         K = [0.1; 0.1]

#         sol = IntervalObservers.solve(sys, K, f_plus, f_minus, x0_plus, x0_minus, tspan; x0=x0)

#         n = sys.n
#         x = get_state(sol, n)
#         xl = get_lower_nonlinear(sol, n)
#         xu = get_upper_nonlinear(sol, n)

#         # Check that true state is contained in the interval (with tolerance for numerical errors)
#         tolerance = 1e-4
#         for t_idx in 1:length(sol.t)
#             @test all(xl[:, t_idx] .<= x[:, t_idx] .+ tolerance)
#             @test all(x[:, t_idx] .<= xu[:, t_idx] .+ tolerance)
#         end
#     end

#     @testset "Interval Width Evolution" begin
#         A = [-2.0 1.0;
#              0.5 -3.0]
#         C = [1.0, 0.0]

#         f_plus = [(t, y) -> 0.1*y, (t, y) -> 0.05*y]
#         f_minus = [(t, y) -> 0.05*y, (t, y) -> 0.02*y]

#         sys = NonLinearSystem(A, C, f_plus, f_minus; check_metzler=false)
#         tspan = (0.0, 3.0)

#         x0_minus = [0.1, 0.2]
#         x0_plus = [0.3, 0.5]

#         K = [0.15; 0.15]  # Stronger feedback

#         sol = IntervalObservers.solve(sys, K, f_plus, f_minus, x0_plus, x0_minus, tspan)

#         n = sys.n
#         xl = get_lower_nonlinear(sol, n)
#         xu = get_upper_nonlinear(sol, n)

#         # Initial interval width
#         initial_width = sum(xu[:, 1] .- xl[:, 1])

#         # Final interval width (should be narrower due to observer convergence)
#         final_width = sum(xu[:, end] .- xl[:, end])

#         @test final_width <= initial_width + 1e-5  # Allow small numerical tolerance
#     end
end

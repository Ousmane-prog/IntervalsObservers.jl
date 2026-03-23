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
        A = [-2.0 -1.0;  
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
        C = [1.0, 0.0, 0.5]  
        
        f_plus = [(t, y) -> 0.1*y, (t, y) -> 0.05*y]
        f_minus = [(t, y) -> 0.05*y, (t, y) -> 0.02*y]
        
        @test_throws DimensionMismatchError NonLinearSystem(A, C, f_plus, f_minus)
    end

    # @testset "Nonlinear Functions Mismatch" begin
    #     A = [-2.0 1.0;
    #          0.5 -3.0]
    #     C = [1.0, 0.0]
        
    #     f_plus = [(t, y) -> 0.1*y, (t, y) -> 0.05*y]  # 2 functions
    #     f_minus = [(t, y) -> 0.05*y]  # Only 1 function
        
    #     @test_throws DimensionMismatchError NonLinearSystem(A, C, f_plus, f_minus)
    # end
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

        # @test sol.t[1] == tspan[1]
        # @test sol.t[end] == tspan[2]
        # @test size(sol.u, 1) == 3  
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

        # @test sol.t[1] == tspan[1]
        # @test sol.t[end] == tspan[2]
    end
end

@testset "Testing LinearSystem design" begin
 
    @testset "LinearSystem Construction" begin
        A = [-2.0 1.0;
            1.0 -3.0]
        
        C = [1.0, 0.0]

        sys = LinearSystem(A, C)

        @test sys.A == A 
        @test sys.C == C
        @test sys.observable == true
        @test sys.is_metzler == true
    end 
    
    # @testset "Non observable System" begin
    #     A = [-1.0 0.0;
    #         0.0 -2.0]
    #     C = [0.0, 1.0]

    #     sys = LinearSystem(A, C)

    #     @test sys.observable == false
        
    # end

    @testset "Dimension Mismatch" begin
        A = [-1.0 0.0;
            0.0 -2.0]
        C = [0.0, 1.0, 2.5]

        @test_throws DimensionMismatchError validate_system_dimensions(A, C)
    end 

    @testset "Testing Interval Observer for Linear Systems" begin
        A = [-3.0  1.0  0.5;
              1.0 -4.0  1.0;
              0.5  1.0 -2.0]

        C = [1.0, 0.0, 0.0]

        x0  = [1.0, 0.5, 0.3]     
        xl0 = [0.6, 0.2, 0.0]      
        xu0 = [1.4, 0.8, 0.6] 

        sys = LinearSystem(A, C)
        
        tspan = (0.0, 3.0)

        K = positive_interval_gain(sys)
        
        sol = IntervalObservers.solve(
            sys, x0, xl0, xu0,
            K, tspan
        )

        plot_state_intervals(sol, sys)
    end
end 


@testset "Testing NonLinearSystem design" begin

    @testset "Testing Interval Observer for Non Linear Systems" begin
        A = [-3.0  1.0  0.5;
              1.0 -4.0  1.0;
              0.5  1.0 -2.0]
        
        α₁ = 0.5
        m1 = 0.1
        α₂ = 0.3
        m2 = 0.1
        m3 = 0.05
        β  = 1.0

        β₁ = α₁ + m1           
        β₂ = α₂ + m2          
        β₃ = m3           

        A = [-β₁ 0.0 0.0;
             α₁ -β₂ 0.0;
             0.0 α₂ -β₃]

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
        tspan = (0.0, 15.0)

        x0_minus = [0.05; 0.1; 0.15]    
        x0_plus  = [0.15; 0.3; 0.45]    


        # K = [0.1; 0.1; 0.2]
        K = positive_interval_gain(sys, desired_poles = [-1.0, -2.0, -3.0])
        @test length(K) == sys.n
        # @test K == [0.0; 0.0; 0.1]

        sol = IntervalObservers.solve(sys, K, f_plus, f_minus, x0_plus, x0_minus, tspan)
        
        @test sol.t[1] == tspan[1]
        @test sol.t[end] ≈ tspan[2]
    end
end 




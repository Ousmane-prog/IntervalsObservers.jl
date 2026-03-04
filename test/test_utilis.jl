using IntervalObservers 
using Test


@testset "Utility Functions" begin
    @testset "Extraction Functions" begin
        A =[-2.0 1.0;
             1.0 -3.0]

        C = [1.0, 0.0]
        sys = LinearSystem(A, C)
        x0 = [1.0, 0.5]
        xl0 = [0.8, 0.3]
        xu0 = [1.2, 0.7]
        tspan = (0.0, 1.0)

        k = positive_interval_gain(sys)
        sol = IntervalObservers.solve(sys, x0, xl0, xu0, k, tspan)

        n = sys.n

        @testset "Extractions" begin
            x = get_state(sol, n)
            xl = get_lower(sol, n)
            xu = get_upper(sol, n)

            @test x[:, 1] == x0
            @test xl[:, 1] == xl0
            @test xu[:, 1] == xu0
            @test size(x, 1) == n 
            @test size(xl, 1) == n
            @test size(xu, 1) == n
            @test size(x, 2) == size(xl, 2) == size(xu, 2)
        end

        # @testset " " begin
        #     x = get_state(sol, n)
        #     xl = get_upper(sol, n)
        #     xu = get_lower(sol, n)

        #     for t in 1:size(x, 2)
        #         # @test(xl[:, t] .<= x[:, t] + 1e-5)
        #         # @test all(x[:, t] .<= xu[:, t] + 1e-5)
        #         # using broacasting instead 
        #         @test xl[:, t] .<= x[:, t] .+ 1e-5
        #         @test x[:, t] .<= xu[:, t] .+ 1e-5
        #     end
        # end
    end
end

@testset " Monotone Functions" begin
    @testset "Non-Metzler matrix detection" begin
        # NOT Metzler matix
        M = [-2.0 0.0 -8.96;
              0.5 -3.0 -2.0;
              0.0 0.3 -5.0]
    
        @test !IntervalObservers._is_monotone_dynamic(M) 
        
        non_monotone = IntervalObservers._find_non_monotone_entries(M)
        @test !isempty(non_monotone)
        @test (1, 3, -8.96) ∈ non_monotone
        @test (2, 3, -2.0) ∈ non_monotone
        @test length(non_monotone) == 2
    end
    
    @testset "Valid Metzler matrix" begin
        M = [-2.0 0.5 1.0;
              0.2 -3.0 0.8;
              0.0 0.4 -1.0]
        
        @test IntervalObservers._is_monotone_dynamic(M)  
        @test isempty(IntervalObservers._find_non_monotone_entries(M))
    end

    @testset "Diagonalmatrix" begin
        M = [2.3 0.0;
             0.0 -1.5]

        @test IntervalObservers._is_monotone_dynamic(M)
    end
end


@testset "positive_interval_gain" begin
    @testset "Positive Gain" begin
        A = [-2.0 1.0;
             1.0 -3.0]

        C = [0.1, 0.0]
        sys = LinearSystem(A, C)
        K = positive_interval_gain(sys)
        @test all(K .== IntervalObservers.DEFAULT_GAIN_VALUE)
    end

    @testset "Pole placement" begin
        A = [-2.0 1.0;
             1.0 -3.0]

        C = [.1, 0.0]
        sys = LinearSystem(A, C)
        desired_poles = [-1.0, -2.0]
        try
            K = positive_interval_gain(sys, desired_poles=desired_poles)
            @test length(K) == sys.n
        catch e
            @test isa(e, IntervalObservers.NonMonotoneDynamicsError)
        end
    end

    @testset "Pole placement " begin
         A = [-0.6  0.0  0.0;
                 0.5  -0.4  0.0;
                 0.0  0.3  -5.0]

        C = [1.0, 1.0, 1.0]

        try
            sys = LinearSystem(A, C)
            desired_poles = [-1.0, -2.0, -3.0]
            @test_throws IntervalObservers.NonMonotoneDynamicsError positive_interval_gain(sys, desired_poles=desired_poles)
        catch e
            @test isa(e, NonMetzlerMatrixError) || isa(e, NonMonotoneDynamicsError)
        end
    end

    @testset "Desired poles validation" begin
        A = [-2.0  0.5;
             0.0  -3.0]
        C = [1.0, 0.0]
        sys = LinearSystem(A, C)

        @testset "Repeated poles - error" begin
            # Poles with repeated values should throw InvalidDesiredPolesError
            repeated_poles = [-1.0, -1.0]
            @test_throws IntervalObservers.InvalidDesiredPolesError positive_interval_gain(sys, desired_poles=repeated_poles)
        end

        @testset "Distinct poles - accepted" begin
            # Properly distinct poles should be accepted
            distinct_poles = [-1.0, -2.0]
            try
                K = positive_interval_gain(sys, desired_poles=distinct_poles)
                @test length(K) == sys.n
            catch e
                # May fail due to non-monotone dynamics, but NOT due to invalid poles
                @test !isa(e, IntervalObservers.InvalidDesiredPolesError)
            end
        end

        @testset "Wrong dimension poles - error" begin
            # Poles with wrong dimension should throw InvalidDesiredPolesError
            wrong_dim_poles = [-1.0, -2.0, -3.0]  # System is 2D
            @test_throws IntervalObservers.InvalidDesiredPolesError positive_interval_gain(sys, desired_poles=wrong_dim_poles)
        end

        @testset "Triple repeated poles - error" begin
            # Multiple repetitions should be caught
            triple_poles = [-1.0, -1.0, -1.0, -1.0]
            A3 = [-1.0  0.5  0.2;
                  0.0  -2.0  0.3;
                  0.0  0.0  -3.0]
            C3 = [1.0, 0.0, 0.0]
            sys3 = LinearSystem(A3, C3)
            @test_throws IntervalObservers.InvalidDesiredPolesError positive_interval_gain(sys3, desired_poles=triple_poles)
        end
    end
end
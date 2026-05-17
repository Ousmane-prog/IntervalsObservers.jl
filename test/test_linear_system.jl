using IntervalObservers
using Test

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

    @testset "Dimension Mismatch" begin
        A = [-1.0 0.0;
              0.0 -2.0]

        C = [0.0, 1.0, 2.5]

        @test_throws DimensionMismatchError validate_system_dimensions(A, C)
    end

    # @testset "Testing Interval Observer for Linear Systems" begin
    #     A = [-3.0  1.0  0.5;
    #           1.0 -4.0  1.0;
    #           0.5  1.0 -2.0]

    #     C = [1.0, 0.0, 0.0]

    #     x0  = [1.0, 0.5, 0.3]
    #     xl0 = [0.6, 0.2, 0.0]
    #     xu0 = [1.4, 0.8, 0.6]

    #     sys = LinearSystem(A, C)

    #     tspan = (0.0, 3.0)

    #     K = positive_interval_gain(sys)

    #     sol = IntervalObservers.solve(sys, K, x0, xl0, xu0, tspan)

    #     @test sol.t[1] == tspan[1]
    #     @test sol.t[end] == tspan[2]

    #     Z = reduce(hcat, sol.u)

    #     n = sys.n

    #     x  = get_state(Z, n)
    #     xl = get_lower(Z, n)
    #     xu = get_upper(Z, n)

    #     @test x[:, 1] ≈ x0
    #     @test xl[:, 1] ≈ xl0
    #     @test xu[:, 1] ≈ xu0

    #     @test size(x, 1) == n
    #     @test size(xl, 1) == n
    #     @test size(xu, 1) == n
    # end
end


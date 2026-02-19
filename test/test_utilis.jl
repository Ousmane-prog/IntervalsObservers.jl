using IntervalObservers
@testset " Test GETER function" begin
    @testset " Extraction consistency" begin
        A = [-2.0 1.0;
              1.0 -3.0]

        C = [1.0, 0.0]
        sys = LinearSystem(A, C)

        x0 = [1.0, 0.5]
        xl0 = [0.8, 0.3]
        xu0 = [1.2, 0.7]
        tspan = (0.0, 1.0)
        K = positive_interval_gain(sys)
        sol = IntervalObservers.solve(
            sys,
            x0,
            xl0,
            xu0,
            K,
            tspan
        )

        n = sys.n
        x = get_state(sol, n)
        xl = get_lower(sol, n)
        xu = get_upper(sol, n)

        @test x[:,1] == x0
        @test xl[:,1] == xl0
        @test xu[:,1] == xu0

    end
end
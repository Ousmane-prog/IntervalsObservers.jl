@testset "Testing NonLinearSystem design" begin

    @testset "Testing Interval Observer for Nonlinear Systems" begin

        α₁ = 0.5
        m1 = 0.1

        α₂ = 0.3
        m2 = 0.1

        m3 = 0.05

        β₁ = α₁ + m1
        β₂ = α₂ + m2
        β₃ = m3

        A = [
            -β₁ 0.0 0.0;
             α₁ -β₂ 0.0;
             0.0  α₂ -β₃
        ]

        C = [0.0, 0.0, 1.0]

        function f₁(t, a, y)
            return a * y / (1.0 + y)
        end

        function f₃(t, c, y)
            return -c * y
        end

        a_max = 0.4
        a_min = 0.1

        c_max = 0.1
        c_min = 0.01

        f_plus = [
            (t, y) -> f₁(t, a_max, y),
            (t, y) -> 0.0,
            (t, y) -> f₃(t, c_max, y),
        ]

        f_minus = [
            (t, y) -> 0.5 * f₁(t, a_min, y),
            (t, y) -> 0.0,
            (t, y) -> 0.5 * f₃(t, c_min, y),
        ]

        f_true = [
            (t, y) -> f₁(t, 0.25, y),
            (t, y) -> 0.0,
            (t, y) -> f₃(t, 0.05, y),
        ]

        sys = NonLinearSystem(A, C, f_plus, f_minus)

        tspan = (0.0, 5.0)

        x0_minus = [0.05, 0.1, 0.15]
        x0_plus  = [0.15, 0.3, 0.45]
        x0       = [0.1, 0.2, 0.3]

        @test sys.observable == true

        K = positive_interval_gain(
            sys;
            desired_poles = [-1.0, -2.0, -3.0],
        )

        @test length(K) == sys.n

        sol = IntervalObservers.solve(
            sys,
            K,
            f_plus,
            f_minus,
            x0_plus,
            x0_minus,
            tspan;
            x0 = x0,
            f_true = f_true,
        )

        @test sol isa IntervalObservers.IntervalObserverSolution

        @test sol.t[1] == tspan[1]
        @test sol.t[end] == tspan[2]

        Z = reduce(hcat, sol.u)

        n = sys.n

        @test size(Z, 1) == 3n

        x  = get_state(Z, n)
        xl = get_lower(Z, n)
        xu = get_upper(Z, n)

        @test size(x, 1) == n
        @test size(xl, 1) == n
        @test size(xu, 1) == n

        @test x[:, 1] ≈ x0
    end

    @testset "Testing σ solve method" begin

        α₁ = 0.5
        m1 = 0.1

        α₂ = 0.3
        m2 = 0.1

        m3 = 0.05

        β₁ = α₁ + m1
        β₂ = α₂ + m2
        β₃ = m3

        A = [
            -β₁ 0.0 0.0;
             α₁ -β₂ 0.0;
             0.0  α₂ -β₃
        ]

        C = [0.0, 0.0, 1.0]

        f_plus = [
            (t, y) -> 0.4 * y / (1.0 + y),
            (t, y) -> 0.0,
            (t, y) -> -0.01 * y,
        ]

        f_minus = [
            (t, y) -> 0.1 * y / (1.0 + y),
            (t, y) -> 0.0,
            (t, y) -> -0.1 * y,
        ]

        sys = NonLinearSystem(A, C, f_plus, f_minus)

        x0_minus = [0.05, 0.1, 0.15]
        x0_plus  = [0.15, 0.3, 0.45]
        x0       = [0.1, 0.2, 0.3]

        tspan = (0.0, 2.0)

        sol = IntervalObservers.solve(
            sys,
            2.0,
            f_plus,
            f_minus,
            x0_plus,
            x0_minus,
            tspan;
            x0 = x0,
        )

        @test sol isa IntervalObservers.IntervalObserverSolution
        @test length(sol.t) > 0

        Z = reduce(hcat, sol.u)

        @test size(Z, 1) == 9
    end

    @testset "Testing _solve_with_sigma" begin

        α₁ = 0.5
        m1 = 0.1

        α₂ = 0.3
        m2 = 0.1

        m3 = 0.05

        β₁ = α₁ + m1
        β₂ = α₂ + m2
        β₃ = m3

        A = [
            -β₁ 0.0 0.0;
             α₁ -β₂ 0.0;
             0.0  α₂ -β₃
        ]

        C = [0.0, 0.0, 1.0]

        function f₁(t, a, y)
            return a * y / (1.0 + y)
        end

        function f₃(t, c, y)
            return -c * y
        end

        f_plus = [
            (t, y) -> f₁(t, 0.4, y),
            (t, y) -> 0.0,
            (t, y) -> f₃(t, 0.01, y),
        ]

        f_minus = [
            (t, y) -> f₁(t, 0.1, y),
            (t, y) -> 0.0,
            (t, y) -> f₃(t, 0.1, y),
        ]

        f_true = [
            (t, y) -> f₁(t, 0.25, y),
            (t, y) -> 0.0,
            (t, y) -> f₃(t, 0.05, y),
        ]

        sys = NonLinearSystem(A, C, f_plus, f_minus)

        x0_minus = [0.05, 0.1, 0.15]
        x0_plus  = [0.15, 0.3, 0.45]
        x0       = [0.1, 0.2, 0.3]

        tspan = (0.0, 2.0)
        saveat = collect(range(tspan[1], tspan[2]; length = 11))

        sol = IntervalObservers._solve_with_sigma(
            sys,
            2.0,
            f_plus,
            f_minus,
            x0_plus,
            x0_minus,
            tspan;
            x0 = x0,
            f_true = f_true,
            saveat = saveat,
        )

        @test sol isa IntervalObservers.IntervalObserverSolution
        @test sol.t == saveat
        @test sol.label == "σ = 2.0, poles = [-2.0, -4.0, -8.0]"

        Z = reduce(hcat, sol.u)
        n = sys.n

        @test size(Z, 1) == 3n
        @test get_state(Z, n)[:, 1] ≈ x0
    end

    @testset "Testing λ interval solve method" begin

        α₁ = 0.5
        m1 = 0.1

        α₂ = 0.3
        m2 = 0.1

        m3 = 0.05

        β₁ = α₁ + m1
        β₂ = α₂ + m2
        β₃ = m3

        A = [
            -β₁ 0.0 0.0;
             α₁ -β₂ 0.0;
             0.0  α₂ -β₃
        ]

        C = [0.0, 0.0, 1.0]

        f_plus = [
            (t, y) -> 0.4 * y / (1.0 + y),
            (t, y) -> 0.0,
            (t, y) -> -0.01 * y,
        ]

        f_minus = [
            (t, y) -> 0.1 * y / (1.0 + y),
            (t, y) -> 0.0,
            (t, y) -> -0.1 * y,
        ]

        sys = NonLinearSystem(A, C, f_plus, f_minus)

        x0_minus = [0.05, 0.1, 0.15]
        x0_plus  = [0.15, 0.3, 0.45]
        x0       = [0.1, 0.2, 0.3]

        tspan = (0.0, 2.0)

        sol = IntervalObservers.solve(
            sys,
            (-5.0, -1.0),
            f_plus,
            f_minus,
            x0_plus,
            x0_minus,
            tspan;
            x0 = x0,
        )

        @test sol isa IntervalObservers.IntervalObserverSolution
        @test sol.label == "intersection"

        Z = reduce(hcat, sol.u)

        @test size(Z, 1) == 9
    end
end
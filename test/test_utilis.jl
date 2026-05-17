using IntervalObservers
using Test
using LinearAlgebra

@testset "Utility Functions" begin
    @testset "Extraction Functions" begin
        A = [-2.0 1.0;
              1.0 -3.0]

        C = [1.0, 0.0]
        sys = LinearSystem(A, C)

        x0  = [1.0, 0.5]
        xl0 = [0.8, 0.3]
        xu0 = [1.2, 0.7]
        tspan = (0.0, 1.0)

        K = positive_interval_gain(sys)

        sol = IntervalObservers.solve(sys, K, x0, xl0, xu0, tspan)
        n = sys.n
        Z = hcat(sol.u...)

        x  = get_state(Z, n)
        xl = get_lower(Z, n)
        xu = get_upper(Z, n)

        # @test x[:, 1] ≈ x0
        @test isapprox(x[:, 1], x0; atol = 1e-8)
        @test isapprox(xl[:, 1], xl0; atol = 1e-8)
        @test isapprox(xu[:, 1], xu0; atol = 1e-8)

        @test size(x, 1) == n
        @test size(xl, 1) == n
        @test size(xu, 1) == n
        @test size(x, 2) == size(xl, 2) == size(xu, 2)
    end
end

@testset "Monotone Functions" begin
    @testset "Non-Metzler matrix detection" begin
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

    @testset "Diagonal matrix" begin
        M = [2.3 0.0;
             0.0 -1.5]

        @test IntervalObservers._is_monotone_dynamic(M)
    end
end

@testset "Interval Transformation" begin
    P = [2.0 -1.0;
        -3.0  4.0]

    x_minus = [1.0, 2.0]
    x_plus  = [3.0, 5.0]

    z_minus, z_plus = IntervalObservers.transform_interval(P, x_minus, x_plus)

    @test z_minus ≈ [-3.0, 5.0]
    @test z_plus  ≈ [4.0, 17.0]

    x0 = [2.0, 3.0]
    z_minus2, z_plus2, z0 = IntervalObservers.compute_bounds_in_new_basis(
        P,
        x_minus,
        x_plus,
        x0,
    )

    # @test z_minus2 ≈ z_minus
    # @test z_plus2  ≈ z_plus
    # @test z0 ≈ P * x0
    @test isapprox(z_minus2, z_minus; atol = 1e-8)
    @test isapprox(z_plus2, z_plus; atol = 1e-8)
    @test isapprox(z0, P * x0; atol = 1e-8)
end

@testset "Pole Generation and Companion/Vandermonde Utilities" begin
    @testset "generate_poles_from_sigma" begin
        poles = IntervalObservers.generate_poles_from_sigma(2.0, 3)

        @test poles == [-2.0, -4.0, -8.0]
        @test_throws AssertionError IntervalObservers.generate_poles_from_sigma(1.0, 3)
        @test_throws AssertionError IntervalObservers.generate_poles_from_sigma(-2.0, 3)
    end

    @testset "create_collection and generate_poles" begin
        λ_vals = (-10.0, -1.0)

        collection = IntervalObservers.create_collection(λ_vals, 3)
        @test collection ≈ [-10.0, -5.5, -1.0]

        pole_sets = IntervalObservers.generate_poles(λ_vals, 3; δ = 0.5)

        @test length(pole_sets) == 3
        @test pole_sets[1] ≈ [-10.0, -5.0, -2.5]
        @test pole_sets[2] ≈ [-5.5, -2.75, -1.375]
        @test pole_sets[3] ≈ [-1.0, -0.5, -0.25]
    end

    @testset "vandermonde_matrix" begin
        λ = [-2.0, -4.0, -8.0]

        V = IntervalObservers.vandermonde_matrix(λ)

        @test size(V) == (3, 3)
        @test V[:, 1] ≈ ones(3)
        @test V[:, 2] ≈ λ
        @test V[:, 3] ≈ λ .^ 2
        @test abs(det(V)) > 1e-10
    end

    @testset "companion_change_of_basis and sigma_change_of_basis" begin
        A = [-0.6 0.0 0.0;
              0.5 -0.4 0.0;
              0.0 0.3 -0.05]

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

        P = IntervalObservers.companion_change_of_basis(sys)

        @test size(P) == (3, 3)
        @test abs(det(P)) > 1e-10

        M, M_inv, poles = IntervalObservers.sigma_change_of_basis(sys, 2.0)

        @test poles ≈ [-2.0, -4.0, -8.0]
        @test size(M) == (3, 3)
        @test size(M_inv) == (3, 3)
        @test M * M_inv ≈ I(3) atol = 1e-8
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

        C = [0.1, 0.0]
        sys = LinearSystem(A, C)

        desired_poles = [-1.0, -2.0]

        try
            K = positive_interval_gain(sys, desired_poles = desired_poles)
            @test length(K) == sys.n
        catch e
            @test isa(e, IntervalObservers.NonMonotoneDynamicsError)
        end
    end

    @testset "Desired poles validation" begin
        A = [-2.0 0.5;
              0.0 -3.0]

        C = [1.0, 0.0]
        sys = LinearSystem(A, C)

        @testset "Repeated poles - error" begin
            repeated_poles = [-1.0, -1.0]

            @test_throws IntervalObservers.InvalidDesiredPolesError positive_interval_gain(
                sys;
                desired_poles = repeated_poles,
            )
        end

        @testset "Distinct poles - accepted" begin
            distinct_poles = [-1.0, -2.0]

            try
                K = positive_interval_gain(sys; desired_poles = distinct_poles)
                @test length(K) == sys.n
            catch e
                @test !isa(e, IntervalObservers.InvalidDesiredPolesError)
            end
        end

        @testset "Wrong dimension poles - error" begin
            wrong_dim_poles = [-1.0, -2.0, -3.0]

            @test_throws IntervalObservers.InvalidDesiredPolesError positive_interval_gain(
                sys;
                desired_poles = wrong_dim_poles,
            )
        end

        @testset "Triple repeated poles - error" begin
            triple_poles = [-1.0, -1.0, -1.0, -1.0]

            A3 = [-1.0 0.5 0.2;
                   0.0 -2.0 0.3;
                   0.0 0.0 -3.0]

            C3 = [1.0, 0.0, 0.0]
            sys3 = LinearSystem(A3, C3)

            @test_throws IntervalObservers.InvalidDesiredPolesError positive_interval_gain(
                sys3;
                desired_poles = triple_poles,
            )
        end
    end
end

@testset "Nonlinear solve methods" begin
    A = [-0.6 0.0 0.0;
          0.5 -0.4 0.0;
          0.0 0.3 -0.05]

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

    tspan = (0.0, 0.5)
    saveat = collect(range(tspan[1], tspan[2]; length = 11))

    @testset "solve with K" begin
        K = positive_interval_gain(sys; desired_poles = [-1.0, -2.0, -3.0])

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
            saveat = saveat,
        )

        @test sol isa IntervalObservers.IntervalObserverSolution
        @test length(sol.t) == length(saveat)
        @test length(sol.u) == length(saveat)

        Z = hcat(sol.u...)
        @test size(Z, 1) == 9
        @test size(Z, 2) == length(saveat)
    end

    @testset "solve with scalar sigma" begin
        sol = IntervalObservers.solve(
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
        @test length(sol.t) == length(saveat)

        Z = hcat(sol.u...)
        @test size(Z, 1) == 9
    end

    @testset "solve with one-value tuple sigma" begin
        sol = IntervalObservers.solve(
            sys,
            (2.0,),
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
        @test length(sol.t) == length(saveat)

        Z = hcat(sol.u...)
        @test size(Z, 1) == 9
    end

    @testset "solve with lambda interval returns intersection" begin
        sol = IntervalObservers.solve(
            sys,
            (-4.0, -1.0),
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
        @test sol.label == "intersection"
        @test length(sol.t) == length(saveat)

        Z = hcat(sol.u...)
        @test size(Z, 1) == 9
    end
end
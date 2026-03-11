using IntervalObservers
using Test
using LinearAlgebra

@testset "IntervalObservers.jl" begin
    # Write your tests here.
    include("test_linear_system.jl")
    include("test_utilis.jl")
    include("test_nonlinear_system.jl")
end

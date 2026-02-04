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

end 
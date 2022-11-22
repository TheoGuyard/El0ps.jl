@testset "Data" begin
    @testset "Regression" begin
        x, A, y = synthetic_data_regression(5, 50, 100, 0.1, 10., normalize=true)
        @test all(norm(ai) ≈ 1. for ai in eachcol(A))
        @test norm(A * x)^2 / norm(y - A * x)^2 ≈ 10.
        @test_throws AssertionError synthetic_data_regression(200, 50, 100, 0.1, 10.)
        @test_throws AssertionError synthetic_data_regression(50, 50, 100, -1., 10.)
        @test_throws AssertionError synthetic_data_regression(50, 50, 100, 0.1, -1.)
    end
    @testset "Classification" begin
        x, A, y = synthetic_data_classification(5, 50, 100, 0.1, 10., normalize=true)
        @test all(norm(ai) ≈ 1. for ai in eachcol(A))
        @test all(abs.(y) .== 1.)
        @test_throws AssertionError synthetic_data_classification(200, 50, 100, 0.1, 10.)
        @test_throws AssertionError synthetic_data_classification(50, 50, 100, -1., 10.)
        @test_throws AssertionError synthetic_data_classification(50, 50, 100, 0.1, -1.)
    end
end
@testset "Bounding" begin

    struct NewBoundingSolver <: El0ps.AbstractBoundingSolver end

    problem = El0ps.Problem(El0ps.LeastSquares(randn(3)), El0ps.Bigm(1.0), randn(3, 5), 0.1)
    bounding_solver = NewBoundingSolver()

    @test_throws ErrorException El0ps.bounding_type(bounding_solver)
    @test_throws ErrorException El0ps.bound!(
        bounding_solver,
        problem,
        nothing,
        nothing,
        nothing,
    )
end

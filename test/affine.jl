@testset "Affine interval matrix in one interval" begin
    A0 = [1 0; 0 1.]
    A1 = [0 1; 1 0.]
    λ = 0 .. 1
    P = AffineIntervalMatrix1(A0, A1, λ)

    @test size(P) == (2, 2)
    @test P[1, 1] == interval(1)
    @test P[1, 2] == interval(0, 1)
    P[1, 1] = 2.0
    @test P[1, 1] == interval(2)
end

@testset "Affine interval matrix in several intervals" begin
    A0 = [1 0; 0 1.]
    A1 = [0 1; 1 0.]; A2 = copy(A1)
    λ1 = 0 .. 1; λ2 = copy(λ1)
    P = AffineIntervalMatrix(A0, [A1, A2], [λ1, λ2])

    @test size(P) == (2, 2)
    @test P[1, 1] == interval(1)
    @test P[1, 2] == interval(0, 2)
    P[1, 1] = 5.0
    @test P[1, 1] == interval(5)
end

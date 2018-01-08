using IntervalMatrices, Base.Test

@testset "Interval arithmetic" begin
    a = -1.5 ± 0.5
    b = -1..1
    @test a + b == -3.0..0.0
    @test a * b == -2.0..2.0
    @test a * b + a == -4.0..1.0
    @test a * (b + 1) == -4.0..0.0
end

@testset "Interval matrix arithmetic" begin
    a = 1.0..1.3; b = 2.0..3.4; c = -0.5 ± 0.5; d = 0 ± 0.1
    A = IntervalMatrix([a b; c d])

    # can add, TODO: add test
    B = A + A
    @test B isa IntervalMatrix

    # can multiply, TODO: add test
    B = A * A
    @test B isa IntervalMatrix
end

@testset "Interval matrix exponential" begin
    m = IntervalMatrix([-1.1..0.9 -4.1.. -3.9; 3.9..4.1 -1.1..0.9])
    me = IntervalMatrixExponential(m, 2)
    @test me isa IntervalMatrixExponential
end

using IntervalMatrices, Test, LinearAlgebra

@testset "Interval arithmetic" begin
    a = -1.5 ± 0.5
    b = -1..1
    @test a + b == -3.0..0.0
    @test a * b == -2.0..2.0
    @test a * b + a == -4.0..1.0
    @test a * (b + 1) == -4.0..0.0
end

@testset "Interval matrix arithmetic" begin
    a = 1.0..1.3; b = 2.0..3.4; c = -0.5 ± 0.5; d = 0.0 ± 0.1
    A = IntervalMatrix([a b; c d])

    # can add, TODO: add test
    B = A + A
    @test B isa IntervalMatrix

    # can multiply, TODO: add test
    B = A * A
    @test B isa IntervalMatrix

    # multiply interval and interval matrix
    x = 0.0..2.0
    @test x * A == A * x ==
          IntervalMatrix([0.0..2.6 0.0..6.8; -2.0..0.0 -0.2..0.2])
end

@testset "Interval matrix methods" begin
    m = IntervalMatrix([-1.1..0.9 -4.1.. -3.9; 3.8..4.2 0.0..0.9])
    @test opnorm(m) == opnorm(m, Inf) ≈ 5.2
    @test opnorm(m, 1) ≈ 5.3
    l = inf(m)
    r = sup(m)
end

@testset "Interval matrix exponential" begin
    m = IntervalMatrix([-1.1..0.9 -4.1.. -3.9; 3.9..4.1 -1.1..0.9])
    a_over = expm_overapproximation(m, 1.0, 4)
    a_under = expm_underapproximation(m, 1.0, 4)
end

@testset "Interval matrix split" begin
    m = IntervalMatrix([-1.1..0.9 -4.1.. -3.9; 3.9..4.1 -1.1..0.9])
    c, s = split(m)
    @test c ≈ [-0.1 -4.; 4. -0.1]
    @test s ≈ [1. 0.1; 0.1 1.]
end

@testset "Interval matrix membership" begin
    m = IntervalMatrix([-1.1..0.9 -4.1.. -3.9; 3.9..4.1 -1.1..0.9])
    a1 = [0. -4.;  4. 0.]
    a2 = [0. -3.; 4. 0.]
    @test a1 ∈ m
    @test a2 ∉ m
end

@testset "Interval matrix rand" begin
    m = IntervalMatrix([-1.1..0.9 -4.1.. -3.9; 3.9..4.1 -1.1..0.9])
    for i in 1:3
        @test rand(m) ∈ m
    end
end

@testset "Interval matrix rand" begin
    m = IntervalMatrix([-1.1..0.9 -4.1.. -3.9; 3.9..4.1 -1.1..0.9])
    f = IntervalMatrices.correction_hull(m, 1e-3, 5)
end

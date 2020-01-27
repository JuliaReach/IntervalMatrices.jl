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
    a = 1.0..1.3; b = 2.0..3.5; c = -0.5 ± 0.5; d = 0.0 ± 0.1
    a₊ = 2.0..2.6; b₊ = 4.0..7.0; c₊ = -2.0..0.0; d₊ = -0.2..0.2
    a₋ = -0.3..0.3; b₋ = -1.5..1.5; c₋ = -1.0..1.0; d₋ = -0.2..0.2
    A = IntervalMatrix([a b; c d])

    B = A + A
    @test B isa IntervalMatrix && B == IntervalMatrix([a₊ b₊; c₊ d₊])

    B = A - A
    @test B isa IntervalMatrix && B == IntervalMatrix([a₋ b₋; c₋ d₋])

    # can multiply, TODO: add test
    B = A * A
    @test B isa IntervalMatrix

    # multiply interval and interval matrix
    x = 0.0..2.0
    @test x * A == A * x ==
          IntervalMatrix([0.0..2.6 0.0..7.0; -2.0..0.0 -0.2..0.2])
    # multiply scalar and interval matrix
    x = 1.0
    for A2 in [x * A, A * x]
        @test A2 == A && typeof(A2) == typeof(A)
    end
end

@testset "Interval matrix methods" begin
    m = IntervalMatrix([-1.1..0.9 -4.1.. -3.9; 3.8..4.2 0.0..0.9])
    @test opnorm(m) == opnorm(m, Inf) ≈ 5.2
    @test opnorm(m, 1) ≈ 5.3
    l = inf(m)
    r = sup(m)
    c = mid(m)
    d = diam(m)
    m2 = copy(m)
    @test m2 isa IntervalMatrix && m.mat == m2.mat
    @test l == inf.(m) && r == sup.(m) && c == mid.(m)
    @test d ≈ r - l
end

@testset "Interval matrix exponential" begin
    m = IntervalMatrix([-1.1..0.9 -4.1.. -3.9; 3.9..4.1 -1.1..0.9])
    a_over = expm_overapproximation(m, 1.0, 4)
    a_under = expm_underapproximation(m, 1.0, 4)
    @test a_over isa IntervalMatrix && a_under isa IntervalMatrix
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

@testset "Interval matrix correction terms" begin
    m = IntervalMatrix([-1.1..0.9 -4.1.. -3.9; 3.9..4.1 -1.1..0.9])
    f = IntervalMatrices.correction_hull(m, 1e-3, 5)
    f2 = IntervalMatrices.input_correction(m, 1e-3, 5)
end

@testset "Interval matrix square" begin
    m = IntervalMatrix([-1.1..0.9 -4.1.. -3.9; 3.9..4.1 -1.1..0.9])
    a = m * m
    b = square(m)
    @test all(inf(a) .<= inf(b)) && all(sup(b) .>= sup(a))
end

@testset "Binary decomposition" begin
    bd = IntervalMatrices.binary_decomposition
    DF = IntervalMatrices.DoubledFactor
    OF = IntervalMatrices.OneFactor{Int}()
    @test Dict(1 => [OF]) == bd(1)
    @test Dict(2 => [DF(1)], 1 => [OF]) == bd(2)
    @test Dict(3 => [DF(1), OF], 1 => [OF]) == bd(3)
    @test Dict(9 => [DF(4), OF], 4 => [DF(2)], 2 => [DF(1)], 1 => [OF]) == bd(9)
    @test Dict(65 => [DF(32), OF], 32 => [DF(16)], 16 => [DF(8)], 8 => [DF(4)],
               4 => [DF(2)], 2 => [DF(1)], 1 => [OF]) == bd(65)
end

@testset "Matrix power" begin
    # flat intervals
    B = IntervalMatrix([2.0±0 0.0±0; 0.0±0 3.0±0])
    C = [2.0 0.0; 0.0 3.0]
    @test all(mid(B^k) == C^k for k in 1:20)

    # proper intervals
    D = IntervalMatrix([2.0±0.1 0.0±0.1; 0.0±0.1 3.0±0.1])
    D^4
end

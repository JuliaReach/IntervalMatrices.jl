using IntervalMatrices, Test, LinearAlgebra

using IntervalMatrices: _truncated_exponential_series,
                        horner,
                        scale_and_square,
                        correction_hull,
                        input_correction

@testset "Interval arithmetic" begin
    a = -1.5 ± 0.5
    b = -1..1
    @test a + b == -3.0..0.0
    @test a * b == -2.0..2.0
    @test a * b + a == -4.0..1.0
    @test a * (b + 1) == -4.0..0.0
end

@testset "Interval matrix construction" begin
    m1 = IntervalMatrix([-1.1..0.9 -4.1.. -3.9; 3.8..4.2 0.0..0.9])
    m2 = IntervalMatrix{Float64}(undef, 2, 2)
    @test m2 isa IntervalMatrix{Float64} && size(m2) == (2, 2)
    m3 = similar(m1)
    @test m3 isa IntervalMatrix{Float64} && size(m3) == size(m1)
    m = [1.0 2.0; 3.0 4.0]
    mint = IntervalMatrix([Interval(1) Interval(2); Interval(3) Interval(4)])
    @test IntervalMatrix(m) == mint

    A = [1 2; 3 4]
    B = [1 2; 4 5]

    @test IntervalMatrix(A, B) == IntervalMatrix([1..1 2..2; 3..4 4..5])

    @test A ± B == IntervalMatrix([0..2 0..4;-1..7 -1..9])
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

    B = A * A
    @test B isa IntervalMatrix

    # multiply interval and interval matrix
    x = 0.0..2.0
    @test x * A == A * x ==
          IntervalMatrix([0.0..2.6 0.0..7.0; -2.0..0.0 -0.2..0.2])
    # multiply scalar and interval matrix
    x = 1.0
    for A2 in [x * A, A * x, A / x]
        @test A2 == A && typeof(A2) == typeof(A)
    end

    # arithmetic closure using interval matrices and non-interval matrices
    Ainf = inf(A)
    @test A + Ainf isa IntervalMatrix
    @test Ainf + A isa IntervalMatrix
    @test A - Ainf isa IntervalMatrix
    @test Ainf - A isa IntervalMatrix
    @test A * Ainf isa IntervalMatrix
    @test Ainf * A isa IntervalMatrix
end

@testset "Interval matrix methods" begin
    m = IntervalMatrix([-1.1..0.9 -4.1.. -3.9; 3.8..4.2 0.0..0.9])
    @test opnorm(m) == opnorm(m, Inf) ≈ 5.2
    @test opnorm(m, 1) ≈ 5.3
    l = inf(m)
    r = sup(m)
    c = mid(m)
    d = diam(m)
    rad = radius(m)
    m2 = copy(m)
    @test m2 isa IntervalMatrix && m.mat == m2.mat
    @test l == inf.(m) && r == sup.(m) && c == mid.(m)
    @test d ≈ r - l
    @test rad ≈ d/2
    sm = scale(m, 2.0)
    @test sm ==  2.0 .* m
    @test sm ≠ m
    scale!(m, 2.0) # in-place
    @test sm == m
    m3 = IntervalMatrix([-2.0..2.0 -2.0..0.0; 0.0..2.0 -1.0..1.0])
    m4 = IntervalMatrix([-1.0..1.0 -1.0..1.0; -1.0..1.0 -2.0..2.0])
    @test m3 ∩ m4 == IntervalMatrix([-1.0..1.0 -1.0..0.0; 0.0..1.0 -1.0..1.0])
    @test m3 ∪ m4 == IntervalMatrix([-2.0..2.0 -2.0..1.0; -1.0..2.0 -2.0..2.0])
    @test diam_norm(m3) ≈ 6.0 # default diameter p-norm is Inf
    @test diam_norm(m3, 1) ≈ 6.0
end

@testset "Interval matrix exponential" begin

    @test quadratic_expansion(-3..3, 1.0, 2.0) == Interval(-0.125, 21)

    m = IntervalMatrix([-1.1..0.9 -4.1.. -3.9; 3.9..4.1 -1.1..0.9])

    for i in 0:4
        _truncated_exponential_series(m, 1.0, i)
    end

    overapp1 = exp_overapproximation(m, 1.0, 4)
    overapp2 = horner(m, 10)
    overapp3 = scale_and_square(m, 5, 1.0, 4)
    underapp = exp_underapproximation(m, 1.0, 4)

    @test underapp isa IntervalMatrix
    for overapp in [overapp1, overapp2, overapp3]
        @test overapp isa IntervalMatrix
    end
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

@testset "Interval matrix inclusion" begin
    m1 = IntervalMatrix([-1.1..0.9 -4.1.. -3.9; 3.9..4.1 -1.1..0.9])
    m2 = IntervalMatrix([-1.2..1.0 -4.1.. -3.9; 3.9..4.2 -1.2..0.9])
    @test m1 ⊆ m2 && !(m2 ⊆ m1)

end

@testset "Interval matrix rand" begin
    m = IntervalMatrix([-1.1..0.9 -4.1.. -3.9; 3.9..4.1 -1.1..0.9])
    for i in 1:3
        @test rand(m) ∈ m
    end
end

@testset "Interval matrix correction terms" begin
    m = IntervalMatrix([-1.1..0.9 -4.1.. -3.9; 3.9..4.1 -1.1..0.9])
    f = correction_hull(m, 1e-3, 5)
    f2 = input_correction(m, 1e-3, 5)
    f = correction_hull(mid(m), 1e-3, 5)
end

@testset "Interval matrix square" begin
    m = IntervalMatrix([-1.1..0.9 -4.1.. -3.9; 3.9..4.1 -1.1..0.9])
    a = m * m
    b = square(m)
    @test all(inf(a) .<= inf(b)) && all(sup(b) .>= sup(a))
end

@testset "Interval matrix power" begin
    m = IntervalMatrix([2.0..2.0 2.0..3.0; 0.0..0.0 -1.0..1.0])
    pow = IntervalMatrixPower(m)
    increment!(pow)  # power of two
    increment!(pow, algorithm="multiply")
    increment!(pow)  # power of two
    increment!(pow, algorithm="power")
    increment!(pow, algorithm="decompose_binary")
    increment!(pow, algorithm="intersect")
end

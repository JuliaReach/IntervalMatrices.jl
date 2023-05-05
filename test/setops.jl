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
    @test rad ≈ d / 2

    sm = scale(m, 2.0)
    @test sm == 2.0 .* m
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

@testset "Interval matrix membership" begin
    m = IntervalMatrix([-1.1..0.9 -4.1.. -3.9; 3.9..4.1 -1.1..0.9])
    a1 = [0.0 -4; 4 0]
    a2 = [0.0 -3; 4 0]
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

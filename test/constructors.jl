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

@testset "Interval matrix midpoint_radius" begin
    m = IntervalMatrix([-1.1..0.9 -4.1.. -3.9; 3.9..4.1 -1.1..0.9])
    c, s = midpoint_radius(m)
    @test c ≈ [-0.1 -4.; 4. -0.1]
    @test s ≈ [1. 0.1; 0.1 1.]
end

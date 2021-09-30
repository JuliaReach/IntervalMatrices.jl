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

@testset "Matrix multiplication" begin

    # test default settings
    @test get_multiplication_mode() == Dict(:multiplication => :fast)

    A = IntervalMatrix([2..4 -2..1; -1..2 2..4])
    set_multiplication_mode(:slow)
    @test A * A == IntervalMatrix([0..18 -16..8; -8..16 0..18])
    @test A * mid.(A) == IntervalMatrix([5..12.5 -8..2; -2..8 5..12.5])
    @test mid.(A) * A == IntervalMatrix([5..12.5 -8..2; -2..8 5..12.5])

    # set_multiplication_mode(:rank1)
    # @test A * A == [0..18 -16..8; -8..16 0..18]

    set_multiplication_mode(:fast)
    @test A * A == IntervalMatrix([-2..19.5 -16..10; -10..16 -2..19.5])
    @test A * mid.(A) == IntervalMatrix([5..12.5 -8..2; -2..8 5..12.5])
    @test mid.(A) * A == IntervalMatrix([5..12.5 -8..2; -2..8 5..12.5])
end

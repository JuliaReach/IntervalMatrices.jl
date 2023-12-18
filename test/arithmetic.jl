@testset "Interval arithmetic" begin
    a = interval(-2, -1)
    b = interval(-1, 1)
    @test a + b == interval(-3, 0)
    @test a * b == interval(-2, 2)
    @test a * b + a == interval(-4, 1)
    @test a * (b + 1) == interval(-4, 0)
end

@testset "Interval matrix arithmetic" begin
    a = interval(1, 1.3)
    b = interval(2, 3.5)
    c = interval(-1, 0)
    d = interval(-0.1, 0.1)
    a₊ = interval(2, 2.6)
    b₊ = interval(4, 7)
    c₊ = interval(-2, 0)
    d₊ = interval(-0.2, 0.2)
    a₋ = interval(-0.3, 0.3)
    b₋ = interval(-1.5, 1.5)
    c₋ = interval(-1, 1)
    d₋ = interval(-0.2, 0.2)
    A = IntervalMatrix([a b; c d])

    B = A + A
    @test B isa IntervalMatrix && B == IntervalMatrix([a₊ b₊; c₊ d₊])

    B = A - A
    @test B isa IntervalMatrix && B[1, 1].lo ≈ a₋.lo && B[1, 1].hi ≈ a₋.hi &&
          B[1, 2] == b₋ && B[2, 1] == c₋ && B[2, 2] == d₋

    # arithmetic with an interval or number
    for x in (interval(0, 2), 2.0)
        @test A + x == x + A == IntervalMatrix([a+x b+x; c+x d+x])
    end

    # multiply interval and interval matrix
    x = interval(0, 2)
    Ax = IntervalMatrix([interval(0, 2.6) interval(0, 7); interval(-2, 0) interval(-0.2, 0.2)])
    @test x * A == A * x == Ax
    # multiply scalar and interval matrix
    x = 1.0
    for A2 in [x * A, A * x, A / x]
        @test A2 == A && typeof(A2) == typeof(A)
    end

    # matrix division
    A2 = IntervalMatrix([1.0 0; 0 2])
    B = A2 \ A
    @test typeof(B) == typeof(A) && B == IntervalMatrix([a b; c/2 d/2])

    # arithmetic closure using interval matrices and non-interval matrices
    Ainf = inf(A)
    @test A + Ainf isa IntervalMatrix
    @test Ainf + A isa IntervalMatrix
    @test A - Ainf isa IntervalMatrix
    @test Ainf - A isa IntervalMatrix
    @test A * Ainf isa IntervalMatrix
    @test Ainf * A isa IntervalMatrix
    @test A \ Ainf isa IntervalMatrix
    @test Ainf \ A isa IntervalMatrix
end

@testset "Matrix multiplication" begin

    # test default settings
    @test get_multiplication_mode() == Dict(:multiplication => :fast)

    A = IntervalMatrix([interval(2, 4) interval(-2, 1); interval(-1, 2) interval(2, 4)])
    set_multiplication_mode(:slow)
    @test A * A == IntervalMatrix([interval(0, 18) interval(-16, 8); interval(-8, 16) interval(0, 18)])
    @test A * mid.(A) == IntervalMatrix([interval(5, 12.5) interval(-8, 2); interval(-2, 8) interval(5, 12.5)])
    @test mid.(A) * A == IntervalMatrix([interval(5, 12.5) interval(-8, 2); interval(-2, 8) interval(5, 12.5)])

    # set_multiplication_mode(:rank1)
    # @test A * A == [interval(0, 18) interval(-16, 8); interval(-8, 16) interval(0, 18)]

    set_multiplication_mode(:fast)
    @test A * A == IntervalMatrix([interval(-2, 19.5) interval(-16, 10); interval(-10, 16) interval(-2, 19.5)])
    @test A * mid.(A) == IntervalMatrix([interval(5, 12.5) interval(-8, 2); interval(-2, 8) interval(5, 12.5)])
    @test mid.(A) * A == IntervalMatrix([interval(5, 12.5) interval(-8, 2); interval(-2, 8) interval(5, 12.5)])
end

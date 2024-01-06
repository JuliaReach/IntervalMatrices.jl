using IntervalMatrices: TaylorOverapproximation,
                        TaylorUnderapproximation,
                        ScaleAndSquare,
                        Horner,
                        horner,
                        scale_and_square,
                        _truncated_exponential_series,
                        _exp_remainder_series

@testset "Interval matrix exponential" begin
    @test quadratic_expansion(interval(-3, 3), 1.0, 2.0) == interval(-0.125, 21)

    M = IntervalMatrix([interval(-1.1, 0.9) interval(-4.1, -3.9);
                        interval(3.9, 4.1) interval(-1.1, 0.9)])

    for i in 0:4
        _truncated_exponential_series(M, 1.0, i)
    end

    overapp1 = exp_overapproximation(M, 1.0, 4)
    overapp2 = horner(M, 10)
    @test_throws ArgumentError horner(M, 0)
    @test_throws ArgumentError horner(M, 3)
    overapp3 = scale_and_square(M, 5, 1.0, 4)
    @test_throws ArgumentError scale_and_square(M, 1, 10.0, 4)
    overapp4 = _exp_remainder_series(M, 1.0, 5)
    underapp1 = exp_underapproximation(M, 1.0, 4)

    @test all(x -> isa(x, IntervalMatrix),
              [overapp1, overapp2, overapp3, overapp4, underapp1])

    @test exp(M) == exp(M; alg=ScaleAndSquare(5, 4)) # default
    @test overapp1 == exp(M; alg=TaylorOverapproximation(4))
    @test overapp2 == exp(M; alg=Horner(10))
    @test overapp3 == exp(M; alg=ScaleAndSquare(5, 4))
    @test underapp1 == exp(M; alg=TaylorUnderapproximation(4))
end

@testset "Interval matrix correction terms" begin
    m = IntervalMatrix([interval(-1.1, 0.9) interval(-4.1, -3.9);
                        interval(3.9, 4.1) interval(-1.1, 0.9)])
    f = correction_hull(m, 1e-3, 5)
    f2 = input_correction(m, 1e-3, 5)
    f = correction_hull(mid(m), 1e-3, 5)
end

@testset "Interval matrix square" begin
    m = IntervalMatrix([interval(-1.1, 0.9) interval(-4.1, -3.9);
                        interval(3.9, 4.1) interval(-1.1, 0.9)])

    a = m * m
    b = square(m)
    @test b ⊆ a
end

@testset "Interval matrix power" begin
    m = IntervalMatrix([interval(2, 2) interval(2, 3); interval(0, 0) interval(-1, 1)])
    pow = IntervalMatrixPower(m)

    @test base(pow) === m
    @test get(pow) === pow.Mᵏ
    @test index(pow) === 1

    pow2 = copy(pow)
    @test pow2 isa IntervalMatrixPower && get(pow) == get(pow2)

    pow2 = IntervalMatrixPower(m, 2)
    pow3 = increment(pow)
    @test pow3 isa IntervalMatrix && pow3 == get(pow2)

    increment!(pow)  # next step is a power of two
    @test index(pow) == 2
    @test_throws ArgumentError increment!(pow; algorithm="foo")
    @test index(pow) == 2
    increment!(pow; algorithm="multiply")
    @test index(pow) == 3
    increment!(pow)  # next step is a power of two
    @test index(pow) == 4
    increment!(pow; algorithm="power")
    @test index(pow) == 5
    increment!(pow; algorithm="decompose_binary")
    @test index(pow) == 6
    increment!(pow; algorithm="intersect")
    @test index(pow) == 7
end

@testset "Transmission line model" begin
    A = Matrix(transmission_line())
    expA = exp(A)
    @test opnorm(expA, Inf) < 1e-6

    @static if PkgVersion.Version(IntervalMatrices.IntervalArithmetic) < v"0.22"
        # in newer versions, too large numbers lead to empty interval
        Aint = IntervalMatrix(interval.(A))
        expA_int = exp_overapproximation(Aint, 1.0, 10)
        @test opnorm(expA_int, Inf) > 1e122
    end
end

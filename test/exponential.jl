using IntervalMatrices: TaylorOverapproximation,
                        TaylorUnderapproximation,
                        ScaleAndSquare,
                        Horner

@testset "Interval matrix exponential" begin
    @test quadratic_expansion(-3 .. 3, 1.0, 2.0) == Interval(-0.125, 21)

    M = IntervalMatrix([-1.1 .. 0.9 -4.1 .. -3.9; 3.9 .. 4.1 -1.1 .. 0.9])

    for i in 0:4
        _truncated_exponential_series(M, 1.0, i)
    end

    overapp1 = exp_overapproximation(M, 1.0, 4)
    overapp2 = horner(M, 10)
    overapp3 = scale_and_square(M, 5, 1.0, 4)
    underapp1 = exp_underapproximation(M, 1.0, 4)

    @test all(x -> isa(x, IntervalMatrix), [overapp1, overapp2, overapp3, underapp1])

    @test exp(M) == exp(M; alg=ScaleAndSquare(5, 4)) # default
    @test overapp1 == exp(M; alg=TaylorOverapproximation(4))
    @test overapp2 == exp(M; alg=Horner(10))
    @test overapp3 == exp(M; alg=ScaleAndSquare(5, 4))
    @test underapp1 == exp(M; alg=TaylorUnderapproximation(4))
end

@testset "Interval matrix correction terms" begin
    m = IntervalMatrix([-1.1 .. 0.9 -4.1 .. -3.9; 3.9 .. 4.1 -1.1 .. 0.9])
    f = correction_hull(m, 1e-3, 5)
    f2 = input_correction(m, 1e-3, 5)
    f = correction_hull(mid(m), 1e-3, 5)
end

@testset "Interval matrix square" begin
    m = IntervalMatrix([-1.1 .. 0.9 -4.1 .. -3.9; 3.9 .. 4.1 -1.1 .. 0.9])

    a = m * m
    b = square(m)
    @test b ⊆ a
end

@testset "Interval matrix power" begin
    m = IntervalMatrix([2.0 .. 2.0 2.0 .. 3.0; 0.0 .. 0.0 -1.0 .. 1.0])
    pow = IntervalMatrixPower(m)
    increment!(pow)  # power of two
    increment!(pow; algorithm="multiply")
    increment!(pow)  # power of two
    increment!(pow; algorithm="power")
    increment!(pow; algorithm="decompose_binary")
    increment!(pow; algorithm="intersect")
end

@testset "Transmission line model" begin
    A = Matrix(transmission_line())
    expA = exp(A)
    @test opnorm(expA, Inf) < 1e-6
    Aint = IntervalMatrix(interval.(A))
    expA_int = exp_overapproximation(Aint, 1.0, 10)
    @test opnorm(expA_int, Inf) > 1e122
end

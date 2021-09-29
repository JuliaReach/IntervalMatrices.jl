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

@testset "Transmission line model" begin
    A = transmission_line() |> Matrix
    expA = exp(A)
    @test opnorm(expA, Inf)< 1e-6
    Aint = IntervalMatrix(interval.(A))
    expA_int = exp_overapproximation(Aint, 1.0, 10)
    @test opnorm(expA_int, Inf) > 1e122
end

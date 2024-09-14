# ===============================================
# helper methods for IntervalArithmetic interface
# ===============================================
using IntervalMatrices: intersect_interval, in_interval, issubset_interval

@testset "Interval conversion" begin
    x = convert(Interval{Float64}, 1.0)

    y = convert(Interval{Float64}, x)
    # `===` is not applicable here because it just checks value equivalence
    # for (immutable) `Interval`s
    @test y ⩵ x && y isa Interval{Float64}

    y = convert(Interval{Float32}, x)
    @test y ⩵ x && y isa Interval{Float32}
end

@testset "Interval operation `intersect_interval`" begin
    x = interval(1, 3)
    for y in (interval(2, 4), interval(2.0f0, 4.0f0))
        @test intersect_interval(x, y) ⩵ interval(2, 3)
    end
    @test intersect_interval(x, x) ⩵ x
    for y in (interval(4, 5), interval(4.0f0, 5.0f0))
        @test intersect_interval(x, y) ⩵ emptyinterval(Float64)
    end
end

@testset "Interval operation `in_interval`" begin
    x = interval(1, 3)
    for e in (1, 2.0, 3)
        @test in_interval(e, x)
    end
    for e in (0, 0.999, 3.001, 4)
        @test !in_interval(e, x)
    end
end

@testset "Interval operation `issubset_interval`" begin
    x = interval(1, 3)
    for y in (interval(0, 4), x)
        @test issubset_interval(x, y)
    end
    for y in (interval(2, 4), interval(0, 2), interval(2))
        @test !issubset_interval(x, y)
    end
end

@testset "Interval from a complex number" begin
    @test interval(1 + 2im) isa Complex{Interval{Float64}}
end

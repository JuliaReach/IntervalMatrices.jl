# ===============================================
# helper methods for IntervalArithmetic interface
# ===============================================

@testset "Interval conversion" begin
    x = convert(Interval{Float64}, 1.0)

    y = convert(Interval{Float64}, x)
    # `===` is not applicable here because it just checks value equivalence
    # for (immutable) `Interval`s
    @test y ⩵ x && y isa Interval{Float64}

    y = convert(Interval{Float32}, x)
    @test y ⩵ x && y isa Interval{Float32}
end

@testset "Interval operation `±`" begin
    @test 1 ± 2 ⩵ interval(-1, 3)
end

@testset "Interval operation `intersect`" begin
    x = interval(1, 3)
    for y in (interval(2, 4), interval(2f0, 4f0))
        @test x ∩ y ⩵ interval(2, 3)
    end
    @test x ∩ x ⩵ x
    for y in (interval(4, 5), interval(4f0, 5f0))
        @test x ∩ y ⩵ emptyinterval(Float64)
    end
end

@testset "Interval operation `in`" begin
    x = interval(1, 3)
    for e in (1, 2.0, 3)
        @test e ∈ x
    end
    for e in (0, 0.999, 3.001, 4)
        @test e ∉ x
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

@testset "Equality for `Interval` matrix" begin
    M = [interval(1) interval(2); interval(3) interval(4)]
    @test M == [1 2; 3 4]
    M = [interval(1, 2) interval(2, 3); interval(3, 4) interval(4, 5)]
    @test M == M
end

@testset "Interval from a complex number" begin
    @test interval(1+2im) isa Complex{Interval{Float64}}
end

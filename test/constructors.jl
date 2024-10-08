@testset "Interval matrix construction" begin
    m1 = IntervalMatrix([interval(-1.1, 0.9) interval(-4.1, -3.9);
                         interval(3.8, 4.2) interval(0, 0.9)])
    m2 = IntervalMatrix{Float64}(undef, 2, 2)
    @test m2 isa IntervalMatrix{Float64} && size(m2) == (2, 2)
    m3 = similar(m1)
    @test m3 isa IntervalMatrix{Float64} && size(m3) == size(m1)
    m = [1.0 2; 3 4]
    mint = IntervalMatrix([interval(1) interval(2); interval(3) interval(4)])
    @test IntervalMatrix(m) == mint

    A = [1 2; 3 4]
    B = [1 2; 4 5]

    @test IntervalMatrix(A, B) ==
          IntervalMatrix([interval(1, 1) interval(2, 2); interval(3, 4) interval(4, 5)])

    @test A ± B == IntervalMatrix([interval(0, 2) interval(0, 4); interval(-1, 7) interval(-1, 9)])
end

@testset "Interval matrix midpoint_radius" begin
    m = IntervalMatrix([interval(-1.1, 0.9) interval(-4.1, -3.9);
                        interval(3.9, 4.1) interval(-1.1, 0.9)])
    c, s = midpoint_radius(m)
    @test c ≈ [-0.1 -4; 4 -0.1]
    @test s ≈ [1 0.1; 0.1 1]
end

@testset "Complex interval matrices" begin
    m1 = IntervalMatrix([interval(1, 2)+im * interval(3, 4) 1; 2 3])
    @test m1 isa IntervalMatrix{Float64}
    @test eltype(m1) == Complex{Interval{Float64}}

    rp = real(m1)
    ip = imag(m1)
    @test rp isa IntervalMatrix{Float64}
    @test eltype(rp) == Interval{Float64}

    @test ip isa IntervalMatrix{Float64}
    @test eltype(ip) == Interval{Float64}

    m2 = IntervalMatrix([1+im 2+im; 3+im 4+im])
    @test m2 isa IntervalMatrix{Float64}
    @test eltype(m2) == Complex{Interval{Float64}}
end

@testset "Special matrices" begin
    A = rand(IntervalMatrix)
    @test A isa IntervalMatrix

    Aₛ = Symmetric(A)
    @test Aₛ isa IntervalMatrix
    @test Aₛ.mat isa Symmetric
    @test Matrix(Aₛ) isa Matrix

    @static if VERSION >= v"1.3"
        Aₕ = UpperHessenberg(A)
        @test Aₕ isa IntervalMatrix
        @test Aₕ.mat isa UpperHessenberg
        @test Matrix(Aₕ) isa Matrix
    end
end

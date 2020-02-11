using IntervalMatrices, Test, LinearAlgebra, SparseArrays

using IntervalMatrices: _truncated_exponential_series, scale_and_square

@testset "Interval arithmetic" begin
    a = -1.5 ± 0.5
    b = -1..1
    @test a + b == -3.0..0.0
    @test a * b == -2.0..2.0
    @test a * b + a == -4.0..1.0
    @test a * (b + 1) == -4.0..0.0
end

@testset "Interval matrix construction" begin
    m1 = IntervalMatrix([-1.1..0.9 -4.1.. -3.9; 3.8..4.2 0.0..0.9])
    m2 = IntervalMatrix{Float64}(undef, 2, 2)
    @test m2 isa IntervalMatrix{Float64} && size(m2) == (2, 2)
    m3 = similar(m1)
    @test m3 isa IntervalMatrix{Float64} && size(m3) == size(m1)
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

    B = A * A
    @test B isa IntervalMatrix

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

@testset "Interval matrix methods" begin
    m = IntervalMatrix([-1.1..0.9 -4.1.. -3.9; 3.8..4.2 0.0..0.9])
    @test opnorm(m) == opnorm(m, Inf) ≈ 5.2
    @test opnorm(m, 1) ≈ 5.3
    l = inf(m)
    r = sup(m)
    c = mid(m)
    d = diam(m)
    m2 = copy(m)
    @test m2 isa IntervalMatrix && m.mat == m2.mat
    @test l == inf.(m) && r == sup.(m) && c == mid.(m)
    @test d ≈ r - l
    sm = scale(m, 2.0)
    @test sm ==  2.0 .* m
    @test sm ≠ m
    scale!(m, 2.0) # in-place
    @test sm == m
    m3 = IntervalMatrix([-2.0..2.0 -2.0..0.0; 0.0..2.0 -1.0..1.0])
    m4 = IntervalMatrix([-1.0..1.0 -1.0..1.0; -1.0..1.0 -2.0..2.0])
    @test m3 ∩ m4 == IntervalMatrix([-1.0..1.0 -1.0..0.0; 0.0..1.0 -1.0..1.0])
end

@testset "Interval matrix exponential" begin
    m = IntervalMatrix([-1.1..0.9 -4.1.. -3.9; 3.9..4.1 -1.1..0.9])

    for i in 0:4
        _truncated_exponential_series(m, 1.0, i)
    end

    overapp1 = exp_overapproximation(m, 1.0, 4)
    overapp2 = scale_and_square(m, 5, 1.0, 4)
    underapp = exp_underapproximation(m, 1.0, 4)

    @test underapp isa IntervalMatrix
    for overapp in [overapp1, overapp2]
        @test overapp isa IntervalMatrix
    end
end

@testset "Interval matrix split" begin
    m = IntervalMatrix([-1.1..0.9 -4.1.. -3.9; 3.9..4.1 -1.1..0.9])
    c, s = split(m)
    @test c ≈ [-0.1 -4.; 4. -0.1]
    @test s ≈ [1. 0.1; 0.1 1.]
end

@testset "Interval matrix membership" begin
    m = IntervalMatrix([-1.1..0.9 -4.1.. -3.9; 3.9..4.1 -1.1..0.9])
    a1 = [0. -4.;  4. 0.]
    a2 = [0. -3.; 4. 0.]
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

@testset "Interval matrix correction terms" begin
    m = IntervalMatrix([-1.1..0.9 -4.1.. -3.9; 3.9..4.1 -1.1..0.9])
    f = IntervalMatrices.correction_hull(m, 1e-3, 5)
    f2 = IntervalMatrices.input_correction(m, 1e-3, 5)
end

@testset "Interval matrix square" begin
    m = IntervalMatrix([-1.1..0.9 -4.1.. -3.9; 3.9..4.1 -1.1..0.9])
    a = m * m
    b = square(m)
    @test all(inf(a) .<= inf(b)) && all(sup(b) .>= sup(a))
end

@testset "Interval matrix power" begin
    m = IntervalMatrix([-1.1..0.9 -4.1.. -3.9; 3.9..4.1 -1.1..0.9])
    pow = IntervalMatrixPower(m)
    increment!(pow)  # power of two
    increment!(pow, algorithm="multiply")
    increment!(pow)  # power of two
    increment!(pow, algorithm="power")
    increment!(pow, algorithm="intersect")
end

#=
State matrix `A` from the transmission line model [1]; we refer to that paper
for the interpretation of the physical parameters in this matrix. Here `A` is a
sparse, structured, square matrix of order 2η, where η ≥ 2 is a parameter
(optional, default: 2). For any choice of (positive) parameters, this matrix
is stable in the sense of Hurwitz (i.e. its eigenvalues have strictly negative
real-part). A parameter `scale` (optional, default: 1e-9) to rescale the physical
quantities is applied to `A`.

[1] M. Althoff, B. H. Krogh, and O. Stursberg. Analyzing Reachability of
Linear Dynamic Systems with Parametric Uncertainties. Modeling, Design,
and Simulation of Systems with Uncertainties. 2011.
=#
function tline(; η=2, R₀=1.00, Rd₀=10.0, L₀=1e-10, C₀=1e-13 * 4.00, scale=1e-9)
    A₁₁ = zeros(η, η)
    A₁₂ = Bidiagonal(fill(-1/C₀, η), fill(1/C₀, η-1), :U)
    A₂₁ = Bidiagonal(fill(1/L₀, η), fill(-1/L₀, η-1), :L)
    A₂₂ = Diagonal(vcat(-Rd₀/L₀, fill(-R₀/L₀, η-1)))
    A  = [A₁₁ A₁₂; A₂₁ A₂₂]
    return A .* scale
end

@testset "Transmission line test" begin
    A = tline() |> Matrix |> IntervalMatrix
    power(A, 16) # TODO: test that it doesn't blow-up
    exp_overapproximation(A, 1.0, 10) # TODO: test that it doesn't blow-up
end

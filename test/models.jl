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
function transmission_line(; η=2, R₀=1.00, Rd₀=10.0, L₀=1e-10, C₀=1e-13 * 4.00, scale=1e-9)
    A₁₁ = zeros(η, η)
    A₁₂ = Bidiagonal(fill(-1 / C₀, η), fill(1 / C₀, η - 1), :U)
    A₂₁ = Bidiagonal(fill(1 / L₀, η), fill(-1 / L₀, η - 1), :L)
    A₂₂ = Diagonal(vcat(-Rd₀ / L₀, fill(-R₀ / L₀, η - 1)))
    A = [A₁₁ A₁₂; A₂₁ A₂₂]
    return A .* scale
end

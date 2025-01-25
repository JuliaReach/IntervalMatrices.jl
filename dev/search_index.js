var documenterSearchIndex = {"docs":
[{"location":"about/#About","page":"About","title":"About","text":"","category":"section"},{"location":"about/","page":"About","title":"About","text":"This page contains some general information about this project, and recommendations about contributing.","category":"page"},{"location":"about/","page":"About","title":"About","text":"Pages = [\"about.md\"]","category":"page"},{"location":"about/#Contributing","page":"About","title":"Contributing","text":"","category":"section"},{"location":"about/","page":"About","title":"About","text":"If you like this package, consider contributing! You can send bug reports (or fix them and send your code), add examples to the documentation or propose new features.","category":"page"},{"location":"about/","page":"About","title":"About","text":"Below we detail some of the guidelines that should be followed when contributing to this package. Further information can be found in the JuliaReachDevDocs site.","category":"page"},{"location":"about/#Branches","page":"About","title":"Branches","text":"","category":"section"},{"location":"about/","page":"About","title":"About","text":"Each pull request (PR) should be pushed in a new branch with the name of the author followed by a descriptive name, e.g. mforets/my_feature. If the branch is associated to a previous discussion in one issue, we use the name of the issue for easier lookup, e.g. mforets/7.","category":"page"},{"location":"about/#Unit-testing-and-continuous-integration-(CI)","page":"About","title":"Unit testing and continuous integration (CI)","text":"","category":"section"},{"location":"about/","page":"About","title":"About","text":"This project is synchronized with GitHub Actions such that each PR gets tested before merging (and the build is automatically triggered after each new commit). For the maintainability of this project, it is important to make all unit tests pass.","category":"page"},{"location":"about/","page":"About","title":"About","text":"To run the unit tests locally, you can do:","category":"page"},{"location":"about/","page":"About","title":"About","text":"julia> using Pkg\n\njulia> Pkg.test(\"IntervalMatrices\")","category":"page"},{"location":"about/","page":"About","title":"About","text":"We also advise adding new unit tests when adding new features to ensure long-term support of your contributions.","category":"page"},{"location":"about/#Contributing-to-the-documentation","page":"About","title":"Contributing to the documentation","text":"","category":"section"},{"location":"about/","page":"About","title":"About","text":"This documentation is written in Markdown, and it relies on Documenter.jl to produce the HTML layout. To build the docs, run make.jl:","category":"page"},{"location":"about/","page":"About","title":"About","text":"$ julia --color=yes docs/make.jl","category":"page"},{"location":"about/#Credits","page":"About","title":"Credits","text":"","category":"section"},{"location":"about/","page":"About","title":"About","text":"These persons have contributed to IntervalMatrices.jl (in alphabetic order):","category":"page"},{"location":"about/","page":"About","title":"About","text":"Luca Ferranti\nMarcelo Forets\nChristian Schilling","category":"page"},{"location":"lib/types/#Types","page":"Types","title":"Types","text":"","category":"section"},{"location":"lib/types/","page":"Types","title":"Types","text":"This section describes systems types implemented in IntervalMatrices.jl.","category":"page"},{"location":"lib/types/","page":"Types","title":"Types","text":"Pages = [\"types.md\"]\nDepth = 3","category":"page"},{"location":"lib/types/","page":"Types","title":"Types","text":"CurrentModule = IntervalMatrices","category":"page"},{"location":"lib/types/#Abstract-interval-operators","page":"Types","title":"Abstract interval operators","text":"","category":"section"},{"location":"lib/types/","page":"Types","title":"Types","text":"AbstractIntervalMatrix","category":"page"},{"location":"lib/types/#IntervalMatrices.AbstractIntervalMatrix","page":"Types","title":"IntervalMatrices.AbstractIntervalMatrix","text":"AbstractIntervalMatrix{IT} <: AbstractMatrix{IT}\n\nAbstract supertype for interval matrix types.\n\n\n\n\n\n","category":"type"},{"location":"lib/types/#Interval-matrix","page":"Types","title":"Interval matrix","text":"","category":"section"},{"location":"lib/types/","page":"Types","title":"Types","text":"IntervalMatrix","category":"page"},{"location":"lib/types/#IntervalMatrices.IntervalMatrix","page":"Types","title":"IntervalMatrices.IntervalMatrix","text":"IntervalMatrix{T, IT, MT<:AbstractMatrix{IT}} <: AbstractIntervalMatrix{IT}\n\nAn interval matrix i.e. a matrix whose coefficients are intervals. This type is parameterized in the number field, the interval type, and the matrix type.\n\nFields\n\nmat – matrix whose entries are intervals\n\nExamples\n\njulia> A = IntervalMatrix([-1 .. -0.8 0 .. 0; 0 .. 0 -1 .. -0.8])\n2×2 IntervalMatrix{Float64, Interval{Float64}, Matrix{Interval{Float64}}}:\n [-1.0, -0.7999999]   [0.0, 0.0]\n  [0.0, 0.0]         [-1.0, -0.7999999]\n\nAn interval matrix proportional to the identity matrix can be built using the UniformScaling operator from the standard library LinearAlgebra. For example,\n\njulia> using LinearAlgebra\n\njulia> IntervalMatrix(interval(1)*I, 2)\n2×2 IntervalMatrix{Float64, Interval{Float64}, Matrix{Interval{Float64}}}:\n [1.0, 1.0]  [0.0, 0.0]\n [0.0, 0.0]  [1.0, 1.0]\n\nThe number of columns can be specified as a third argument, creating a rectangular m  n matrix such that only the entries in the main diagonal, (1 1) (2 2)   (k k) are specified, where k = min(m n):\n\njulia> IntervalMatrix(interval(-1, 1)*I, 2, 3)\n2×3 IntervalMatrix{Float64, Interval{Float64}, Matrix{Interval{Float64}}}:\n [-1.0, 1.0]   [0.0, 0.0]  [0.0, 0.0]\n  [0.0, 0.0]  [-1.0, 1.0]  [0.0, 0.0]\n\njulia> IntervalMatrix(interval(-1, 1)*I, 3, 2)\n3×2 IntervalMatrix{Float64, Interval{Float64}, Matrix{Interval{Float64}}}:\n [-1.0, 1.0]   [0.0, 0.0]\n  [0.0, 0.0]  [-1.0, 1.0]\n  [0.0, 0.0]   [0.0, 0.0]\n\nAn uninitialized interval matrix can be constructed using undef:\n\njulia> m = IntervalMatrix{Float64}(undef, 2, 2);\n\njulia> typeof(m)\nIntervalMatrix{Float64, Interval{Float64}, Matrix{Interval{Float64}}}\n\nNote that this constructor implicitly uses a dense matrix, Matrix{Float64}, as the matrix (mat) field in the new interval matrix.\n\n\n\n\n\n","category":"type"},{"location":"lib/types/#Interval-matrix-power-wrapper","page":"Types","title":"Interval-matrix-power wrapper","text":"","category":"section"},{"location":"lib/types/","page":"Types","title":"Types","text":"IntervalMatrixPower","category":"page"},{"location":"lib/types/#IntervalMatrices.IntervalMatrixPower","page":"Types","title":"IntervalMatrices.IntervalMatrixPower","text":"IntervalMatrixPower{T}\n\nA wrapper for the matrix power that can be incremented.\n\nFields\n\nM  – the original matrix\nMᵏ – the current matrix power, i.e., M^k\nk  – the current power index\n\nNotes\n\nThe wrapper should only be accessed using the interface functions. The internal representation (such as the fields) are subject to future changes.\n\nExamples\n\njulia> A = IntervalMatrix([interval(2, 2) interval(2, 3); interval(0, 0) interval(-1, 1)])\n2×2 IntervalMatrix{Float64, Interval{Float64}, Matrix{Interval{Float64}}}:\n [2.0, 2.0]   [2.0, 3.0]\n [0.0, 0.0]  [-1.0, 1.0]\n\njulia> pow = IntervalMatrixPower(A);\n\njulia> increment!(pow)\n2×2 IntervalMatrix{Float64, Interval{Float64}, Matrix{Interval{Float64}}}:\n [4.0, 4.0]  [2.0, 9.0]\n [0.0, 0.0]  [0.0, 1.0]\n\njulia> increment(pow)\n2×2 IntervalMatrix{Float64, Interval{Float64}, Matrix{Interval{Float64}}}:\n [8.0, 8.0]  [-1.0, 21.0]\n [0.0, 0.0]  [-1.0, 1.0]\n\njulia> matrix(pow)\n2×2 IntervalMatrix{Float64, Interval{Float64}, Matrix{Interval{Float64}}}:\n [4.0, 4.0]  [2.0, 9.0]\n [0.0, 0.0]  [0.0, 1.0]\n\njulia> index(pow)\n2\n\njulia> base(pow)\n2×2 IntervalMatrix{Float64, Interval{Float64}, Matrix{Interval{Float64}}}:\n [2.0, 2.0]   [2.0, 3.0]\n [0.0, 0.0]  [-1.0, 1.0]\n\n\n\n\n\n\n","category":"type"},{"location":"lib/types/#Affine-interval-matrix","page":"Types","title":"Affine interval matrix","text":"","category":"section"},{"location":"lib/types/","page":"Types","title":"Types","text":"AffineIntervalMatrix1\nAffineIntervalMatrix","category":"page"},{"location":"lib/types/#IntervalMatrices.AffineIntervalMatrix1","page":"Types","title":"IntervalMatrices.AffineIntervalMatrix1","text":"AffineIntervalMatrix1{T, IT, MT0<:AbstractMatrix{T}, MT1<:AbstractMatrix{T}} <: AbstractIntervalMatrix{IT}\n\nInterval matrix representing the matrix\n\nA₀ + λA₁\n\nwhere A₀ and A₁ are real (or complex) matrices, and λ is an interval.\n\nFields\n\nA0 – matrix\nA1 – matrix\nλ  – interval\n\nExamples\n\nThe matrix I + 1 1 -1 1 * interval(0 1) is:\n\njulia> using LinearAlgebra\n\njulia> P = AffineIntervalMatrix1(Matrix(1.0I, 2, 2), [1 1; -1 1.], interval(0, 1));\n\njulia> P\n2×2 AffineIntervalMatrix1{Float64, Interval{Float64}, Matrix{Float64}, Matrix{Float64}}:\n  [1.0, 2.0]  [0.0, 1.0]\n [-1.0, 0.0]  [1.0, 2.0]\n\n\n\n\n\n","category":"type"},{"location":"lib/types/#IntervalMatrices.AffineIntervalMatrix","page":"Types","title":"IntervalMatrices.AffineIntervalMatrix","text":"AffineIntervalMatrix{T, IT, MT0<:AbstractMatrix{T}, MT<:AbstractMatrix{T}, MTA<:AbstractVector{MT}} <: AbstractIntervalMatrix{IT}\n\nInterval matrix representing the matrix\n\nA₀ + λ₁A₁ + λ₂A₂ +  + λₖAₖ\n\nwhere A₀ and A₁  Aₖ are real (or complex) matrices, and λ₁  λₖ are intervals.\n\nFields\n\nA0 – matrix\nA  – vector of matrices\nλ  – vector of intervals\n\nNotes\n\nThis type is the general case of the AffineIntervalMatrix1, which only contains one matrix proportional to an interval.\n\nExamples\n\nThe affine matrix I + 1 1 -1 1 * interval(0 1) + 0 1 1 0 * interval(2 3) is:\n\njulia> using LinearAlgebra\n\njulia> A0 = Matrix(1.0I, 2, 2);\n\njulia> A1 = [1 1; -1 1.]; A2 = [0 1; 1 0];\n\njulia> λ1 = interval(0, 1); λ2 = interval(2, 3);\n\njulia> P = AffineIntervalMatrix(A0, [A1, A2], [λ1, λ2])\n2×2 AffineIntervalMatrix{Float64, Interval{Float64}, Matrix{Float64}, Matrix{Float64}, Vector{Matrix{Float64}}, Vector{Interval{Float64}}}:\n [1.0, 2.0]  [2.0, 4.0]\n [1.0, 3.0]  [1.0, 2.0]\n\n\n\n\n\n","category":"type"},{"location":"bibliography/#Bibliography","page":"Bibliography","title":"Bibliography","text":"","category":"section"},{"location":"bibliography/","page":"Bibliography","title":"Bibliography","text":"M. Althoff. Reachability analysis and its application to the safety assessment of autonomous cars. Ph.D. Thesis, Technische Universität München (2010).\n\n\n\nM. Althoff, B. H. Krogh and O. Stursberg. Analyzing reachability of linear dynamic systems with parametric uncertainties. Modeling, Design, and Simulation of Systems with Uncertainties, 69–94 (2011).\n\n\n\nM. Althoff, O. Stursberg and M. Buss. Reachability analysis of linear systems with uncertain parameters                   and inputs. In: Conference on Decision and Control (CDC) (IEEE, 2007); pp. 726–732.\n\n\n\nM. Althoff, O. Stursberg and M. Buss. Reachability analysis of nonlinear systems with uncertain parameters                   using conservative linearization. In: Conference on Decision and Control (CDC) (IEEE, 2008); pp. 4042–4048.\n\n\n\nA. Goldsztejn and A. Neumaier. On the Exponentiation of Interval Matrices. Reliable Computing 20, 53–72 (2014).\n\n\n\nO. Kosheleva, V. Kreinovich, G. Mayer and H. T. Nguyen. Computing the cube of an interval matrix is NP-Hard. In: Symposium on Applied Computing (SAC), edited by H. Haddad, L. M. Liebrock, A. Omicini and R. L. Wainwright (ACM, 2005); pp. 1449–1453.\n\n\n\nM. Liou. A novel method of evaluating transient response. Proceedings of the IEEE 54, 20–23 (1966).\n\n\n\nS. M. Rump. Verification methods: Rigorous results using floating-point arithmetic. Acta Numerica 19, 287–449 (2010).\n\n\n\nS. M. Rump. Fast and parallel interval arithmetic. BIT Numerical Mathematics 39, 534–554 (1999).\n\n\n\n","category":"page"},{"location":"#IntervalMatrices.jl","page":"Home","title":"IntervalMatrices.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"IntervalMatrices is a Julia package to work with matrices that have interval coefficients.","category":"page"},{"location":"#Features","page":"Home","title":"Features","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Here is a quick summary of the available functionality. See the section Library Outline below for details.","category":"page"},{"location":"","page":"Home","title":"Home","text":"An IntervalMatrix type that wraps a two-dimensional array whose components are intervals.\nArithmetic between scalars, intervals and interval matrices.\nQuadratic expansions using single-use expressions (SUE).\nInterval matrix exponential: underapproximation and overapproximation routines.\nUtility functions such as: operator norm, random sampling, splitting and containment check.","category":"page"},{"location":"","page":"Home","title":"Home","text":"An application of interval matrices is to find the set of states reachable by a dynamical system whose coefficients are uncertain. The library ReachabilityAnalysis.jl implements algorithms that use interval matrices [AKS11, Lio66].","category":"page"},{"location":"#Installing","page":"Home","title":"Installing","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Depending on your needs, choose an appropriate command from the following list and enter it in Julia's REPL.","category":"page"},{"location":"","page":"Home","title":"Home","text":"To install the latest release version:","category":"page"},{"location":"","page":"Home","title":"Home","text":"pkg> add IntervalMatrices","category":"page"},{"location":"","page":"Home","title":"Home","text":"To install the latest development version:","category":"page"},{"location":"","page":"Home","title":"Home","text":"pkg> add IntervalMatrices#master","category":"page"},{"location":"","page":"Home","title":"Home","text":"To clone the package for development:","category":"page"},{"location":"","page":"Home","title":"Home","text":"pkg> dev IntervalMatrices","category":"page"},{"location":"#Quickstart","page":"Home","title":"Quickstart","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"An interval matrix is a matrix whose coefficients are intervals. For instance,","category":"page"},{"location":"","page":"Home","title":"Home","text":"julia> using IntervalMatrices\n\njulia> A = IntervalMatrix([interval(0, 1) interval(1, 2); interval(2, 3) interval(-4, -2)])\n2×2 IntervalMatrix{Float64, Interval{Float64}, Matrix{Interval{Float64}}}:\n [0.0, 1.0]   [1.0, 2.0]\n [2.0, 3.0]  [-4.0, -2.0]","category":"page"},{"location":"","page":"Home","title":"Home","text":"defines an interval matrix A. The type of each coefficient in A is an interval, e.g. its coefficient in position (1 1) is the interval 0 1 over double-precision floating-point numbers:","category":"page"},{"location":"","page":"Home","title":"Home","text":"julia> A[1, 1]\n[0.0, 1.0]\n\njulia> typeof(A[1, 1])\nInterval{Float64}","category":"page"},{"location":"","page":"Home","title":"Home","text":"This library uses the interval arithmetic package IntervalArithmetic.jl to deal with interval computations. For instance, one can compute a multiple of A,","category":"page"},{"location":"","page":"Home","title":"Home","text":"julia> 2A\n2×2 IntervalMatrix{Float64, Interval{Float64}, Matrix{Interval{Float64}}}:\n [0.0, 2.0]   [2.0, 4.0]\n [4.0, 6.0]  [-8.0, -4.0]","category":"page"},{"location":"","page":"Home","title":"Home","text":"Or an interval multiple of A,","category":"page"},{"location":"","page":"Home","title":"Home","text":"julia> interval(-1, 1) * A\n2×2 IntervalMatrix{Float64, Interval{Float64}, Matrix{Interval{Float64}}}:\n [-1.0, 1.0]  [-2.0, 2.0]\n [-3.0, 3.0]  [-4.0, 4.0]","category":"page"},{"location":"","page":"Home","title":"Home","text":"Or compute the square of A,","category":"page"},{"location":"","page":"Home","title":"Home","text":"julia> square(A)\n2×2 IntervalMatrix{Float64, Interval{Float64}, Matrix{Interval{Float64}}}:\n   [2.0, 7.0]   [-8.0, -1.0]\n [-12.0, -2.0]   [6.0, 22.0]","category":"page"},{"location":"","page":"Home","title":"Home","text":"In these cases, the rules of interval arithmetic are used; see the Wikipedia page on interval arithmetic for the relevant definitions and algebraic rules that apply.","category":"page"},{"location":"","page":"Home","title":"Home","text":"However, the straightforward application of the rules of interval arithmetic does not always give the exact result: in general it only gives an overapproximation [ASB08, KKMN05]. To illustrate, suppose that we are interested in the quadratic term At + frac12A^2 t^2, which corresponds to the Taylor-series expansion at order two of e^At - I. Then, at t = 10,","category":"page"},{"location":"","page":"Home","title":"Home","text":"julia> A + 1/2 * A^2\n2×2 IntervalMatrix{Float64, Interval{Float64}, Matrix{Interval{Float64}}}:\n  [1.0, 4.50001]  [-3.0, 2.0]\n [-4.0, 2.50001]  [-1.0, 9.0]","category":"page"},{"location":"","page":"Home","title":"Home","text":"However, that result is not tight. The computation can be performed exactly via single-use expressions implemented in this library:","category":"page"},{"location":"","page":"Home","title":"Home","text":"julia> quadratic_expansion(A, 1.0, 0.5)\n2×2 IntervalMatrix{Float64, Interval{Float64}, Matrix{Interval{Float64}}}:\n  [1.0, 4.50001]  [-2.0, 1.0]\n [-3.0, 1.50001]   [1.0, 7.0]","category":"page"},{"location":"","page":"Home","title":"Home","text":"We now obtain an interval matrix that is strictly included in the one obtained from the naive multiplication.","category":"page"},{"location":"","page":"Home","title":"Home","text":"An overapproximation and an underapproximation method at a given order for e^At, where A is an interval matrix, are also available. See the Methods section for details.","category":"page"},{"location":"#Library-Outline","page":"Home","title":"Library Outline","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Explore the types and methods defined in this library by following the links below, or use the search bar in the left to look for a specific keyword in the documentation.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Pages = [\n    \"lib/types.md\",\n    \"lib/methods.md\"\n]\nDepth = 2","category":"page"},{"location":"lib/methods/#Methods","page":"Methods","title":"Methods","text":"","category":"section"},{"location":"lib/methods/","page":"Methods","title":"Methods","text":"This section describes systems methods implemented in IntervalMatrices.jl.","category":"page"},{"location":"lib/methods/","page":"Methods","title":"Methods","text":"Pages = [\"methods.md\"]\nDepth = 3","category":"page"},{"location":"lib/methods/","page":"Methods","title":"Methods","text":"CurrentModule = IntervalMatrices","category":"page"},{"location":"lib/methods/#Common-functions","page":"Methods","title":"Common functions","text":"","category":"section"},{"location":"lib/methods/","page":"Methods","title":"Methods","text":"inf\nsup\nmid\ndiam\nradius\nmidpoint_radius\nrand\nsample\n∈\n±\n⊆\n∩\n∪\nhull","category":"page"},{"location":"lib/methods/#IntervalArithmetic.inf","page":"Methods","title":"IntervalArithmetic.inf","text":"inf(A::IntervalMatrix{T}) where {T}\n\nReturn the infimum of an interval matrix A, which corresponds to taking the element-wise infimum of A.\n\nInput\n\nA – interval matrix\n\nOutput\n\nA scalar matrix whose coefficients are the infima of each element in A.\n\n\n\n\n\n","category":"function"},{"location":"lib/methods/#IntervalArithmetic.sup","page":"Methods","title":"IntervalArithmetic.sup","text":"sup(A::IntervalMatrix{T}) where {T}\n\nReturn the supremum of an interval matrix A, which corresponds to taking the element-wise supremum of A.\n\nInput\n\nA – interval matrix\n\nOutput\n\nA scalar matrix whose coefficients are the suprema of each element in A.\n\n\n\n\n\n","category":"function"},{"location":"lib/methods/#IntervalArithmetic.mid","page":"Methods","title":"IntervalArithmetic.mid","text":"mid(A::IntervalMatrix{T}) where {T}\n\nReturn the midpoint of an interval matrix A, which corresponds to taking the element-wise midpoint of A.\n\nInput\n\nA – interval matrix\n\nOutput\n\nA scalar matrix whose coefficients are the midpoints of each element in A.\n\n\n\n\n\n","category":"function"},{"location":"lib/methods/#IntervalArithmetic.diam","page":"Methods","title":"IntervalArithmetic.diam","text":"diam(A::IntervalMatrix{T}) where {T}\n\nReturn a matrix whose entries describe the diameters of the intervals.\n\nInput\n\nA – interval matrix\n\nOutput\n\nA matrix B of the same shape as A such that B[i, j] == diam(A[i, j]) for each i and j.\n\n\n\n\n\n","category":"function"},{"location":"lib/methods/#IntervalArithmetic.radius","page":"Methods","title":"IntervalArithmetic.radius","text":"radius(A::IntervalMatrix{T}) where {T}\n\nReturn the radius of an interval matrix A, which corresponds to taking the element-wise radius of A.\n\nInput\n\nA – interval matrix\n\nOutput\n\nA scalar matrix whose coefficients are the radii of each element in A.\n\n\n\n\n\n","category":"function"},{"location":"lib/methods/#IntervalArithmetic.midpoint_radius","page":"Methods","title":"IntervalArithmetic.midpoint_radius","text":"midpoint_radius(A::IntervalMatrix{T}) where {T}\n\nSplit an interval matrix A into two scalar matrices C and S such that A = C + -S S.\n\nInput\n\nA – interval matrix\n\nOutput\n\nA pair (C, S) such that the entries of C are the central points and the entries of S are the (nonnegative) radii of the intervals in A.\n\n\n\n\n\n","category":"function"},{"location":"lib/methods/#Base.rand","page":"Methods","title":"Base.rand","text":"rand(::Type{IntervalMatrix}, m::Int=2, [n]::Int=m;\n     N=Float64, rng::AbstractRNG=GLOBAL_RNG)\n\nReturn a random interval matrix of the given size and numeric type.\n\nInput\n\nIntervalMatrix – type, used for dispatch\nm              – (optional, default: 2) number of rows\nn              – (optional, default: m) number of columns\nrng            – (optional, default: GLOBAL_RNG) random-number generator\n\nOutput\n\nAn interval matrix of size m  n whose coefficients are normally-distributed intervals of type N with mean 0 and standard deviation 1.\n\nNotes\n\nIf this function is called with only one argument, it creates a square matrix, because the number of columns defaults to the number of rows.\n\n\n\n\n\n","category":"function"},{"location":"lib/methods/#IntervalMatrices.sample","page":"Methods","title":"IntervalMatrices.sample","text":"sample(A::IntervalMatrix{T}; rng::AbstractRNG=GLOBAL_RNG) where {T}\n\nReturn a sample of the given random interval matrix.\n\nInput\n\nA   – interval matrix\nm   – (optional, default: 2) number of rows\nn   – (optional, default: 2) number of columns\nrng – (optional, default: GLOBAL_RNG) random-number generator\n\nOutput\n\nAn interval matrix of size m  n whose coefficients are normally-distributed intervals of type N with mean 0 and standard deviation 1.\n\n\n\n\n\n","category":"function"},{"location":"lib/methods/#Base.:∈","page":"Methods","title":"Base.:∈","text":"∈(M::AbstractMatrix, A::AbstractIntervalMatrix)\n\nCheck whether a concrete matrix is an instance of an interval matrix.\n\nInput\n\nM – concrete matrix\nA – interval matrix\n\nOutput\n\ntrue iff M is an instance of A\n\nAlgorithm\n\nWe check for each entry in M whether it belongs to the corresponding interval in A.\n\n\n\n\n\n","category":"function"},{"location":"lib/methods/#IntervalArithmetic.:±","page":"Methods","title":"IntervalArithmetic.:±","text":"±(C::MT, S::MT) where {T, MT<:AbstractMatrix{T}}\n\nReturn an interval matrix such that the center and radius of the intervals is given by the matrices C and S respectively.\n\nInput\n\nC – center matrix\nS – radii matrix\n\nOutput\n\nAn interval matrix M such that M[i, j] corresponds to the interval whose center is C[i, j] and whose radius is S[i, j], for each i and j. That is, M = C + -S S.\n\nNotes\n\nThe radii matrix should be nonnegative, i.e. S[i, j] ≥ 0 for each i and j.\n\nExamples\n\njulia> [1 2; 3 4] ± [1 2; 4 5]\n2×2 IntervalMatrix{Float64, Interval{Float64}, Matrix{Interval{Float64}}}:\n  [0.0, 2.0]   [0.0, 4.0]\n [-1.0, 7.0]  [-1.0, 9.0]\n\n\n\n\n\n","category":"function"},{"location":"lib/methods/#Base.:⊆","page":"Methods","title":"Base.:⊆","text":"⊆(A::AbstractIntervalMatrix, B::AbstractIntervalMatrix)\n\nCheck whether an interval matrix is contained in another interval matrix.\n\nInput\n\nA – interval matrix\nB – interval matrix\n\nOutput\n\ntrue iff A[i, j] ⊆ B[i, j] for all i, j.\n\n\n\n\n\n","category":"function"},{"location":"lib/methods/#Base.:∩","page":"Methods","title":"Base.:∩","text":"∩(A::IntervalMatrix, B::IntervalMatrix)\n\nIntersect two interval matrices.\n\nInput\n\nA – interval matrix\nB – interval matrix (of the same shape as A)\n\nOutput\n\nA new matrix C of the same shape as A such that C[i, j] = A[i, j] ∩ B[i, j] for each i and j.\n\n\n\n\n\n","category":"function"},{"location":"lib/methods/#Base.:∪","page":"Methods","title":"Base.:∪","text":"∪(A::IntervalMatrix, B::IntervalMatrix)\n\nFinds the interval union (hull) of two interval matrices. This is equivalent to hull.\n\nInput\n\nA – interval matrix\nB – interval matrix (of the same shape as A)\n\nOutput\n\nA new matrix C of the same shape as A such that C[i, j] = A[i, j] ∪ B[i, j] for each i and j.\n\n\n\n\n\n","category":"function"},{"location":"lib/methods/#IntervalArithmetic.hull","page":"Methods","title":"IntervalArithmetic.hull","text":"hull(A::IntervalMatrix, B::IntervalMatrix)\n\nFinds the interval hull of two interval matrices. This is equivalent to ∪.\n\nInput\n\nA – interval matrix\nB – interval matrix (of the same shape as A)\n\nOutput\n\nA new matrix C of the same shape as A such that C[i, j] = hull(A[i, j], B[i, j]) for each i and j.\n\n\n\n\n\n","category":"function"},{"location":"lib/methods/#Arithmetic","page":"Methods","title":"Arithmetic","text":"","category":"section"},{"location":"lib/methods/","page":"Methods","title":"Methods","text":"square\nscale\nscale!\nset_multiplication_mode","category":"page"},{"location":"lib/methods/#IntervalMatrices.square","page":"Methods","title":"IntervalMatrices.square","text":"square(A::IntervalMatrix)\n\nCompute the square of an interval matrix.\n\nInput\n\nA – interval matrix\n\nOutput\n\nAn interval matrix equivalent to A * A.\n\nAlgorithm\n\nWe follow Kosheleva et al. [KKMN05], Section 6.\n\n\n\n\n\n","category":"function"},{"location":"lib/methods/#IntervalMatrices.scale","page":"Methods","title":"IntervalMatrices.scale","text":"scale(A::IntervalMatrix{T}, α::T) where {T}\n\nReturn a new interval matrix whose entries are scaled by the given factor.\n\nInput\n\nA – interval matrix\nα – scaling factor\n\nOutput\n\nA new matrix B of the same shape as A such that B[i, j] = α*A[i, j] for each i and j.\n\nNotes\n\nSee scale! for the in-place version of this function.\n\n\n\n\n\n","category":"function"},{"location":"lib/methods/#IntervalMatrices.scale!","page":"Methods","title":"IntervalMatrices.scale!","text":"scale!(A::IntervalMatrix{T}, α::T) where {T}\n\nModifies the given interval matrix, scaling its entries by the given factor.\n\nInput\n\nA – interval matrix\nα – scaling factor\n\nOutput\n\nThe matrix A such that for each i and j, the new value of A[i, j] is α*A[i, j].\n\nNotes\n\nThis is the in-place version of scale.\n\n\n\n\n\n","category":"function"},{"location":"lib/methods/#IntervalMatrices.set_multiplication_mode","page":"Methods","title":"IntervalMatrices.set_multiplication_mode","text":"set_multiplication_mode(multype)\n\nSets the algorithm used to perform matrix multiplication with interval matrices.\n\nInput\n\nmultype – symbol describing the algorithm used\n:slow – uses traditional matrix multiplication algorithm.\n:fast – computes an enclosure of the matrix product using the midpoint-radius            notation of the matrix [Rum10].\n\nnote: :fast option no longer supported\n:fast support was removed in IntervalArithmetic v0.22.\n\nNotes\n\nBy default, :slow is used.\nUsing fast is generally significantly faster, but it may return larger intervals, especially if midpoint and radius have the same order of magnitude   (50% overestimate at most) [Rum99].\n\n\n\n\n\n","category":"function"},{"location":"lib/methods/#Matrix-power","page":"Methods","title":"Matrix power","text":"","category":"section"},{"location":"lib/methods/","page":"Methods","title":"Methods","text":"increment!\nincrement\nmatrix\nbase\nindex","category":"page"},{"location":"lib/methods/#IntervalMatrices.increment!","page":"Methods","title":"IntervalMatrices.increment!","text":"increment!(pow::IntervalMatrixPower; [algorithm=default_algorithm])\n\nIncrement a matrix power in-place (i.e., storing the result in pow).\n\nInput\n\npow       – wrapper of a matrix power (modified in this function)\nalgorithm – (optional; default: default_algorithm) algorithm to compute                the matrix power; available options:\n\"multiply\" – fast computation using * from the previous result\n\"power\" – recomputation using ^\n\"decompose_binary\" – decompose k = 2a + b\n\"intersect\" – combination of \"multiply\"/\"power\"/\"decompose_binary\"\n\nOutput\n\nThe next matrix power, reflected in the modified wrapper.\n\nNotes\n\nIndependent of \"algorithm\", if the index is a power of two, we compute the exact result using squaring.\n\n\n\n\n\n","category":"function"},{"location":"lib/methods/#IntervalMatrices.increment","page":"Methods","title":"IntervalMatrices.increment","text":"increment(pow::IntervalMatrixPower; [algorithm=default_algorithm])\n\nIncrement a matrix power without modifying pow.\n\nInput\n\npow – wrapper of a matrix power\nalgorithm – (optional; default: default_algorithm) algorithm to compute                the matrix power; see increment! for available options\n\nOutput\n\nThe next matrix power.\n\n\n\n\n\n","category":"function"},{"location":"lib/methods/#IntervalMatrices.matrix","page":"Methods","title":"IntervalMatrices.matrix","text":"matrix(pow::IntervalMatrixPower)\n\nReturn the matrix represented by a wrapper of a matrix power.\n\nInput\n\npow – wrapper of a matrix power\n\nOutput\n\nThe matrix power represented by the wrapper.\n\n\n\n\n\n","category":"function"},{"location":"lib/methods/#IntervalMatrices.base","page":"Methods","title":"IntervalMatrices.base","text":"base(pow::IntervalMatrixPower)\n\nReturn the original matrix represented by a wrapper of a matrix power.\n\nInput\n\npow – wrapper of a matrix power\n\nOutput\n\nThe matrix M being the basis of the matrix power M^k represented by the wrapper.\n\n\n\n\n\n","category":"function"},{"location":"lib/methods/#IntervalMatrices.index","page":"Methods","title":"IntervalMatrices.index","text":"index(pow::IntervalMatrixPower)\n\nReturn the current index of the wrapper of a matrix power.\n\nInput\n\npow – wrapper of a matrix power\n\nOutput\n\nThe index k of the wrapper representing M^k.\n\n\n\n\n\n","category":"function"},{"location":"lib/methods/#Matrix-exponential","page":"Methods","title":"Matrix exponential","text":"","category":"section"},{"location":"lib/methods/#Algorithms","page":"Methods","title":"Algorithms","text":"","category":"section"},{"location":"lib/methods/","page":"Methods","title":"Methods","text":"IntervalMatrices.Horner\nIntervalMatrices.ScaleAndSquare\nIntervalMatrices.TaylorOverapproximation\nIntervalMatrices.TaylorUnderapproximation","category":"page"},{"location":"lib/methods/#IntervalMatrices.Horner","page":"Methods","title":"IntervalMatrices.Horner","text":"Horner <: AbstractExponentiationMethod\n\nMatrix exponential using Horner's method.\n\nFields\n\nK – number of expansions in the Horner scheme\n\n\n\n\n\n","category":"type"},{"location":"lib/methods/#IntervalMatrices.ScaleAndSquare","page":"Methods","title":"IntervalMatrices.ScaleAndSquare","text":"ScaleAndSquare <: AbstractExponentiationMethod\n\nFields\n\nl – scaling-and-squaring order\np – order of the approximation\n\n\n\n\n\n","category":"type"},{"location":"lib/methods/#IntervalMatrices.TaylorOverapproximation","page":"Methods","title":"IntervalMatrices.TaylorOverapproximation","text":"TaylorOverapproximation <: AbstractExponentiationMethod\n\nMatrix exponential overapproximation using a truncated Taylor series.\n\nFields\n\np – order of the approximation\n\n\n\n\n\n","category":"type"},{"location":"lib/methods/#IntervalMatrices.TaylorUnderapproximation","page":"Methods","title":"IntervalMatrices.TaylorUnderapproximation","text":"TaylorUnderapproximation <: AbstractExponentiationMethod\n\nMatrix exponential underapproximation using a truncated Taylor series.\n\nFields\n\np – order of the approximation\n\n\n\n\n\n","category":"type"},{"location":"lib/methods/#Implementations","page":"Methods","title":"Implementations","text":"","category":"section"},{"location":"lib/methods/","page":"Methods","title":"Methods","text":"exp_overapproximation\nhorner\nscale_and_square\nexp_underapproximation","category":"page"},{"location":"lib/methods/#IntervalMatrices.exp_overapproximation","page":"Methods","title":"IntervalMatrices.exp_overapproximation","text":"exp_overapproximation(A::IntervalMatrix{T}, t, p)\n\nOverapproximation of the exponential of an interval matrix, exp(A*t), using a truncated Taylor series.\n\nInput\n\nA – interval matrix\nt – exponentiation factor\np – order of the approximation\n\nOutput\n\nA matrix enclosure of exp(A*t), i.e. an interval matrix M = (m_{ij}) such that [exp(A*t)]_{ij} ⊆ m_{ij}.\n\nAlgorithm\n\nSee Althoff et al. [ASB07], Theorem 1.\n\n\n\n\n\n","category":"function"},{"location":"lib/methods/#IntervalMatrices.horner","page":"Methods","title":"IntervalMatrices.horner","text":"horner(A::IntervalMatrix{T}, K::Integer; [validate]::Bool=true)\n\nCompute the matrix exponential using the Horner scheme.\n\nInput\n\nA – interval matrix\nK – number of expansions in the Horner scheme\nvalidate – (optional; default: true) option to validate the precondition               of the algorithm\n\nAlgorithm\n\nWe use the algorithm in Goldsztejn and Neumaier [GN14], Section 4.2.\n\n\n\n\n\n","category":"function"},{"location":"lib/methods/#IntervalMatrices.scale_and_square","page":"Methods","title":"IntervalMatrices.scale_and_square","text":"scale_and_square(A::IntervalMatrix{T}, l::Integer, t, p;\n                 [validate]::Bool=true)\n\nCompute the matrix exponential using scaling and squaring.\n\nInput\n\nA – interval matrix\nl – scaling-and-squaring order\nt – non-negative time value\np – order of the approximation\nvalidate – (optional; default: true) option to validate the precondition               of the algorithm\n\nAlgorithm\n\nWe use the algorithm in Goldsztejn and Neumaier [GN14], Section 4.3, which first scales A by factor 2^-l, computes the matrix exponential for the scaled matrix, and then squares the result l times.\n\n    exp(A * 2^-l)^2^l\n\n\n\n\n\n","category":"function"},{"location":"lib/methods/#IntervalMatrices.exp_underapproximation","page":"Methods","title":"IntervalMatrices.exp_underapproximation","text":"exp_underapproximation(A::IntervalMatrix{T}, t, p) where {T}\n\nUnderapproximation of the exponential of an interval matrix, exp(A*t), using a truncated Taylor series expansion.\n\nInput\n\nA – interval matrix\nt – exponentiation factor\np – order of the approximation\n\nOutput\n\nAn underapproximation of exp(A*t), i.e. an interval matrix M = (m_{ij}) such that m_{ij} ⊆ [exp(A*t)]_{ij}.\n\nAlgorithm\n\nSee Althoff et al. [ASB07], Theorem 2.\n\n\n\n\n\n","category":"function"},{"location":"lib/methods/#Finite-expansions","page":"Methods","title":"Finite expansions","text":"","category":"section"},{"location":"lib/methods/","page":"Methods","title":"Methods","text":"quadratic_expansion","category":"page"},{"location":"lib/methods/#IntervalMatrices.quadratic_expansion","page":"Methods","title":"IntervalMatrices.quadratic_expansion","text":"quadratic_expansion(A::IntervalMatrix, α::Real, β::Real)\n\nCompute the quadratic expansion of an interval matrix, αA + βA^2, using interval arithmetic.\n\nInput\n\nA – interval matrix\nα – linear coefficient\nβ – quadratic coefficient\n\nOutput\n\nAn interval matrix that encloses B = αA + βA^2.\n\nAlgorithm\n\nThis a variation of the algorithm in Kosheleva et al. [KKMN05], Section 6. If A = (aᵢⱼ) and B = αA + βA^2 = (bᵢⱼ), the idea is to compute each bᵢⱼ by factoring out repeated expressions (thus the term single-use expressions).\n\nFirst, let i = j. In this case,\n\nbⱼⱼ = βsum_k k  j a_jk a_kj + (α + βa_jj) a_jj\n\nNow consider i  j. Then,\n\nbᵢⱼ = βsum_k k  i k  j a_ik a_kj + (α + βa_ii + βa_jj) a_ij\n\n\n\n\n\n","category":"function"},{"location":"lib/methods/#Correction-terms","page":"Methods","title":"Correction terms","text":"","category":"section"},{"location":"lib/methods/","page":"Methods","title":"Methods","text":"correction_hull\ninput_correction","category":"page"},{"location":"lib/methods/#IntervalMatrices.correction_hull","page":"Methods","title":"IntervalMatrices.correction_hull","text":"correction_hull(A::IntervalMatrix{T}, t, p) where {T}\n\nCompute the correction term for the convex hull of a point and its linear map with an interval matrix in order to contain all trajectories of a linear system.\n\nInput\n\nA – interval matrix\nt – non-negative time value\np – order of the approximation\n\nOutput\n\nAn interval matrix representing the correction term.\n\nAlgorithm\n\nSee Althoff et al. [ASB07], Theorem 3.\n\n\n\n\n\n","category":"function"},{"location":"lib/methods/#IntervalMatrices.input_correction","page":"Methods","title":"IntervalMatrices.input_correction","text":"input_correction(A::IntervalMatrix{T}, t, p) where {T}\n\nCompute the input correction matrix for discretizing an inhomogeneous affine dynamical system with an interval matrix and an input domain not containing the origin.\n\nInput\n\nA – interval matrix\nt – non-negative time value\np – order of the Taylor approximation\n\nOutput\n\nAn interval matrix representing the correction matrix.\n\nAlgorithm\n\nSee Althoff [Alt10], Proposition 3.4.\n\n\n\n\n\n","category":"function"},{"location":"lib/methods/#Norms","page":"Methods","title":"Norms","text":"","category":"section"},{"location":"lib/methods/","page":"Methods","title":"Methods","text":"opnorm\ndiam_norm","category":"page"},{"location":"lib/methods/#LinearAlgebra.opnorm","page":"Methods","title":"LinearAlgebra.opnorm","text":"opnorm(A::IntervalMatrix, p::Real=Inf)\n\nThe matrix norm of an interval matrix.\n\nInput\n\nA – interval matrix\np – (optional, default: Inf) the class of p-norm\n\nNotes\n\nThe matrix p-norm of an interval matrix A is defined as\n\n    A_p = max(textinf(A) textsup(A))_p\n\nwhere max and  are taken elementwise.\n\n\n\n\n\n","category":"function"},{"location":"lib/methods/#IntervalMatrices.diam_norm","page":"Methods","title":"IntervalMatrices.diam_norm","text":"diam_norm(A::IntervalMatrix, p=Inf)\n\nReturn the diameter norm of the interval matrix.\n\nInput\n\nA – interval matrix\np – (optional, default: Inf) the p-norm used; valid options are:        1, 2, Inf\n\nOutput\n\nThe operator norm, in the p-norm, of the scalar matrix obtained by taking the element-wise diam function, where diam(x) := sup(x) - inf(x) for an interval x.\n\nNotes\n\nThis function gives a measure of the width of the interval matrix.\n\n\n\n\n\n","category":"function"}]
}

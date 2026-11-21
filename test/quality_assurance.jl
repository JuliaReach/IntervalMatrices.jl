using IntervalMatrices, Test
import PkgVersion, Aqua, ExplicitImports

@testset "ExplicitImports tests" begin
    @static if PkgVersion.Version(IntervalMatrices.IntervalArithmetic) >= v"0.22"
        ignores = (:GLOBAL_RNG,)
    else
        ignores = (:GLOBAL_RNG, :checksquare, :numtype)
    end
    @test isnothing(ExplicitImports.check_all_explicit_imports_are_public(IntervalMatrices;
                                                                          ignore=ignores))
    @test isnothing(ExplicitImports.check_all_explicit_imports_via_owners(IntervalMatrices))
    @static if VERSION >= v"1.9"
        ignores = (:HermOrSym, :AdjOrTrans)
    else
        ignores = (:HermOrSym, :AdjOrTrans, :AbstractTriangular, :Version)
    end
    @test isnothing(ExplicitImports.check_all_qualified_accesses_are_public(IntervalMatrices;
                                                                            ignore=ignores))
    @test isnothing(ExplicitImports.check_all_qualified_accesses_via_owners(IntervalMatrices))
    ignores = (:IntervalArithmetic, :Interval, :in_interval, :intersect_interval, :interval,
               :isequal_interval, :issubset_interval, :setdisplay)  # due to reexporting IntervalArithmetic
    @test isnothing(ExplicitImports.check_no_implicit_imports(IntervalMatrices; ignore=ignores))
    @test isnothing(ExplicitImports.check_no_self_qualified_accesses(IntervalMatrices))
    @test isnothing(ExplicitImports.check_no_stale_explicit_imports(IntervalMatrices))
end

import Pkg
@static if VERSION >= v"1.10"
    # JET v0.9.0 (earliest supported version) requires Julia v1.10
    Pkg.add("JET")
    import JET

    @testset "JET tests" begin
        # false positives for Base functionality
        JET.test_package(IntervalMatrices; target_modules=(IntervalMatrices,))
    end
end

@testset "Aqua tests" begin
    # PkgVersion is only used in old versions
    @static if VERSION >= v"1.9"
        stale_deps = (ignore=[:PkgVersion],)
        ambiguities = true
    else
        stale_deps = true
        ambiguities = (broken=true,)
    end

    # old versions of IntervalArithmetic did not define `convert`
    @static if PkgVersion.Version(IntervalMatrices.IntervalArithmetic) >= v"0.22"
        piracies = true
    else
        piracies = (broken=true,)
    end

    Aqua.test_all(IntervalMatrices; stale_deps=stale_deps, ambiguities=ambiguities,
                  piracies=piracies)
end

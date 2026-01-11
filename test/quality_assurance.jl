using IntervalMatrices, Test
import PkgVersion, Aqua

import Pkg
@static if VERSION >= v"1.6"  # TODO make explicit test requirement
    Pkg.add("ExplicitImports")
    import ExplicitImports

    @testset "ExplicitImports tests" begin
        ignores = (:GLOBAL_RNG,)
        @test isnothing(ExplicitImports.check_all_explicit_imports_are_public(IntervalMatrices;
                                                                              ignore=ignores))
        @test isnothing(ExplicitImports.check_all_explicit_imports_via_owners(IntervalMatrices))
        ignores = (:HermOrSym, :AdjOrTrans)
        @test isnothing(ExplicitImports.check_all_qualified_accesses_are_public(IntervalMatrices;
                                                                                ignore=ignores))
        @test isnothing(ExplicitImports.check_all_qualified_accesses_via_owners(IntervalMatrices))
        ignores = (:IntervalArithmetic, :Interval, :in_interval, :intersect_interval, :interval,
                   :isequal_interval, :issubset_interval, :setdisplay)  # due to reexporting IntervalArithmetic
        @test isnothing(ExplicitImports.check_no_implicit_imports(IntervalMatrices; ignore=ignores))
        @test isnothing(ExplicitImports.check_no_self_qualified_accesses(IntervalMatrices))
        @test isnothing(ExplicitImports.check_no_stale_explicit_imports(IntervalMatrices))
    end
end

@testset "Aqua tests" begin
    # PkgVersion is only used in old versions
    @static if VERSION >= v"1.9"
        stale_deps = (ignore=[:PkgVersion],)
    else
        stale_deps = true
    end

    # old versions of IntervalArithmetic had an undefined export
    @static if PkgVersion.Version(IntervalMatrices.IntervalArithmetic) >= v"0.17.6"
        undefined_exports = true
    else
        undefined_exports = false
    end

    # old versions of IntervalArithmetic did not define `interval` for `Complex` inputs
    @static if PkgVersion.Version(IntervalMatrices.IntervalArithmetic) >= v"0.21"
        piracies = true
    else
        piracies = (broken=true,)
    end

    Aqua.test_all(IntervalMatrices; stale_deps=stale_deps, undefined_exports=undefined_exports,
                  piracies=piracies)
end

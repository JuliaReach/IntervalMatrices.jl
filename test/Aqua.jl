using IntervalMatrices, Test
import PkgVersion, Aqua

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
                  ambiguities=false, piracies=piracies)

    # do not warn about ambiguities in dependencies
    Aqua.detect_ambiguities(IntervalMatrices)
end

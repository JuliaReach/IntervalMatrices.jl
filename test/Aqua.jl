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

    Aqua.test_all(IntervalMatrices; stale_deps=stale_deps, undefined_exports=undefined_exports,
                  # the ambiguities should be resolved in the future
                  ambiguities=(broken=true,),
                  # the piracies should be resolved in the future
                  piracies=(broken=true,))
end

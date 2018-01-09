# About

This page contains some general information about this project, and recommendations
about contributing.

```@contents
Pages = ["about.md"]
```

## Contributing

If you like this package, consider contributing! You can send bug reports (or fix them
and send your code), add examples to the documentation or propose new features.

Below we detail some of the guidelines that should be followed when contributing
to this package.

### Branches

Each pull request (PR) should be pushed in a new branch with the name of the author
followed by a descriptive name, e.g. `t/mforets/my_feature`. If the branch is
associated to a previous discussion in one issue, we use the name of the issue for easier
lookup, e.g. `t/mforets/7`.

### Unit testing and continuous integration (CI)

This project is synchronized with Travis CI, such that each PR gets tested
before merging (and the build is automatically triggered after each new commit).

For the maintainability of this project, we try to understand and fix the failing
doctests if they exist. We develop in Julia v0.6.0, but for experimentation
we also build on the nightly branch.

To run the unit tests locally, you should do:

```julia
$ julia --color=yes test/runtests.jl
```

### Contributing to the documentation

This documentation is written in Markdown, and it relies on
[Documenter.jl](https://github.com/JuliaDocs/Documenter.jl) to produce the HTML
layout. To build the docs, run `make.jl`:

```julia
$ julia --color=yes docs/make.jl
```

## Running the benchmarks

This package contains a suite of benchmarks that is handled through
[PkgBenchmark](https://github.com/JuliaCI/PkgBenchmark.jl). To run the benchmarks,
execute the following commands in Julia's REPL:

```julia
julia> using IntervalMatrices, PkgBenchmark
julia> benchmarkpkg("IntervalMatrices")
```

To save the results to a custom directory, use:

```julia
julia> benchmarkpkg("IntervalMatrices", resultsdir="/Users/forets/Projects")
```

For further options see the inline help of `benchmarkpkg`.

## Credits

These persons have contributed to `IntervalMatrices.jl` (in alphabetic order):

- [Marcelo Forets](http://marcelo-forets.fr)
- [Alexandre Rocca](http://www-verimag.imag.fr/~rocca/)

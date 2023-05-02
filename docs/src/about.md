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
to this package. Further information can be found in the
[`JuliaReachDevDocs` site](https://juliareach.github.io/JuliaReachDevDocs/latest/).

### Branches

Each pull request (PR) should be pushed in a new branch with the name of the author
followed by a descriptive name, e.g. `mforets/my_feature`. If the branch is
associated to a previous discussion in one issue, we use the name of the issue for easier
lookup, e.g. `mforets/7`.

### Unit testing and continuous integration (CI)

This project is synchronized with GitHub Actions such that each PR gets tested
before merging (and the build is automatically triggered after each new commit).
For the maintainability of this project, it is important to make all unit tests
pass.

To run the unit tests locally, you can do:

```julia
julia> using Pkg

julia> Pkg.test("IntervalMatrices")
```

We also advise adding new unit tests when adding new features to ensure
long-term support of your contributions.

### Contributing to the documentation

This documentation is written in Markdown, and it relies on
[Documenter.jl](https://github.com/JuliaDocs/Documenter.jl) to produce the HTML
layout. To build the docs, run `make.jl`:

```bash
$ julia --color=yes docs/make.jl
```

## Credits

These persons have contributed to `IntervalMatrices.jl` (in alphabetic order):

- [Luca Ferranti](https://lucaferranti.github.io/)
- [Marcelo Forets](http://github.com/mforets)
- [Christian Schilling](https://www.christianschilling.net/)

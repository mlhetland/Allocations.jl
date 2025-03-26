# Reference

```@meta
CurrentModule = Allocations
```

```jldoctest ref
julia> using Allocations
```

**Function naming:** Allocation functions that use a straightforward procedure,
or simply use a solver to enforce some property, are named after that procedure
or property (such as `alloc_rand` or `alloc_ef1`). For published algorithms, the
package uses a naming scheme based on the original publication.

The root of the function name is `alloc_`, followed by a publication code:

- For a single-author paper, the first three letters of the author's last name
  are used;
- for multi-author papers, the first letter of the first four authors are
  concatenated.
- To this, the last two digits of the year are added.

For example, the 2/3-MMS algorithm of Garg, McGlaughlin and Taki (2018) is
implemented by `alloc_gmt18`.

- If the same code applies to multiple publications, they are distinguished by
  a suffix `a`, `b`, etc., after the year digits.

- If a single publication discusses multiple algorithms, a number such as `_1`,
  `_2`, etc., is added. So, for example, the third algorithm described by Biswas
  and Barman (2018) is `alloc_bb18_3`.

Some functions (such as `alloc_hh22_1`) are given generic names as aliases (in
this case, `alloc_half_mms`).

!!! note

    These publication codes are similar to the [authorship trigraphs used in
    some citation styles](https://en.wikipedia.org/wiki/Citation). Specifically,
    they follow the conventions of
    [`alpha.bst`](http://tug.ctan.org/tex-archive/biblio/bibtex/base/alpha.bst),
    except that an "et al." character is not added when there are five or more
    authors.

**Accessors:** Using accessors such as `attr(obj)` rather than `obj.attr` comes
with the cost of claiming names in the global namespace, so they are mainly used
when they provide some value, and not as a general convention. Such added value
may be, for example, that it makes the interface clearer or more consistent, or
that it provides some modest level of encapsulation. In many cases, it may also
be worthwhile to use a function to access elements of an internal collection
(such as `foo(X, i)` rather than `X.foo[i]`).

If, however, an object is essentially a named tuple, attributes will tend to be
used without accessors. This applies to thin wrapper (such as `Conflicts`, for
example), where the inner object is simply accessed directly.

Note that some such objects used accessors in earlier versions. These accessors
are still available, but their use is deprecated. (Run `julia` with
`--depwarn=true` to get the appropriate warnings.)


## [Basic types](@id types)

```@autodocs
Modules = [Allocations]
Pages = ["types.jl", "smallbitset.jl"]
```

## Checks and measures

```@autodocs
Modules = [Allocations]
Pages = ["checks.jl", "measures.jl"]
```

## Allocation algorithms

```@autodocs
Modules = [Allocations]
Pages = ["algorithms.jl"]
```

### Reductions

```@autodocs
Modules = [Allocations]
Pages = ["reductions.jl"]
```

## [MIP-based allocation](@id MIP-based)

```@autodocs
Modules = [Allocations]
Pages = ["mip.jl"]
```

## Instance generation

```@autodocs
Modules = [Allocations]
Pages = ["data.jl"]
```

## Matroids

```@autodocs
Modules = [Allocations]
Pages = ["matroids.jl"]
```

## Configuration

```@autodocs
Modules = [Allocations]
Pages = ["conf.jl"]
```

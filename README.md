# Allocations.jl

The Allocations package deals with the fair allocation of indivisible items
to a set of agents. (For more information about this topic, see, e.g., the
[Wikipedia entry on fair item
allocation](https://en.wikipedia.org/wiki/Fair_item_allocation), or [the survey
by Amanatidis et al.](https://arxiv.org/abs/2208.08782))

You can install Allocations.jl using the repository URL. Just start `julia`,
press `]` and add the package:

```
julia> ]
(@v1.8) pkg> add https://github.com/mlhetland/Allocations.jl
```

This should install both the package and its dependencies. For more
documentation, import the package and consult its docstrings. E.g.:

```
julia> using Allocations
...
julia> ?
help?> Allocations
```

More detailed information may be found in the docstrings of individual
functions:

```
help?> alloc_mms
```

Allocation functions that use a straightforward procedure, or simply use
a solver to enforce some property, are named after that procedure or property
(such as `alloc_rand` or `alloc_ef1`). For published algorithms,
the package uses a naming scheme based on the original publication.[^1]
The root of the function name is `alloc_`, followed by a publication code: For
a single-author paper, the first three letters of the author's last name are
used; for multi-author papers, the first letter of the first four authors are
concatenated. To this, the last two digits of the year are added. 

For example, the 2/3-MMS algorithm of Garg, McGlaughlin and Taki (2018) is
implemented by `alloc_gmt18`. (If the same code applies to multiple
publications, they are distinguished by a suffix `a`, `b`, etc., after the year
digits.) If a single publication discusses multiple algorithms, a number such as
`_1`, `_2`, etc., is added. So, for example, the third algorithm described by
Biswas and Barman (2018) is `alloc_bb18_3`.

Some functions (such as `alloc_hh22_1`) are given generic names as aliases (in
this case, `alloc_half_mms`).

[^1]: The publication codes are similar to the [authorship trigraphs used in some citation styles](https://en.wikipedia.org/wiki/Citation). Specifically, they follow the conventions of [`alpha.bst`](http://tug.ctan.org/tex-archive/biblio/bibtex/base/alpha.bst), except that an "et al." character is not added when there are five or more authors.

using Distributions
using Graphs

# Note: Although `Uniform` represents a uniform distribution on [0,1], the
# implementation (https://github.com/JuliaStats/Distributions.jl/blob/a9b0e3c9)
# actually samples from [0,1), by using rand().


"""
    rand_additive(; n=2:10, m=n->2n:4n, v=(n, m)->Uniform(), rng=default_rng())
    rand_profile(; n=2:10, m=n->2n:4n, v=(n, m)->Uniform(), rng=default_rng())

Generate a random additive valuation profile (an `Additive` object) with the
number of agents and items, agents and values chosen using `rand` with the
given values as the first argument. Here `n` is used directly, while `m` should
be a *function of the number of agents*, which returns an argument for `rand`.
Similarly, `v` should be a function of both the number of agents and items.

The defaults for `n` and `m` are taken from [Hummel and
Hetland](https://arxiv.org/abs/2104.06280), who based them on based on
real-world data from [Caragiannis et al.](https://doi.org/10.1145/3355902), with
`n` at most 10, and the average `m/n` approximately 3.

The distribution for the values can be univariate, in which case it is used
independently for each value. For example, to generate Gaussian valuations, with
the [`Distributions`](https://github.com/JuliaStats/Distributions.jl) package,
use `v=(n, m)->Normal()`.

However, it can also be multivariate, in which case each sample should be a
vector of the same length as the number of items, representing the valuation
function of a single agent. For example, to get a `Dirichlet` distribution,
where an agent's values sum to `1`, you could use `v=(n, m)->Dirichlet(m, 2)`.

It is also possible to use a matrix-variate distribution, such as a matrix
normal distribution, where each sample should then be an `m` by `n` matrix, with
each column representing the valuation function of a single agent.

!!! warning

    Note that the matrix samples should be the *transpose* of the ones used in
    the resulting profile. This is to maintain consistency with the multivariate
    distributions, which produce column vectors.

`rand_profile` is an alias for `rand_additive`.
"""
function rand_additive(; n=2:10, m=n->2n:4n, v=(n, m)->Uniform(), rng=default_rng())
    nn = rand(rng, n)
    mm = rand(rng, m(nn))
    vv = v(nn, mm)
    X = rand(rng, nn, mm)
    rand!(rng, vv, X')
    return Additive(X)
end

const rand_profile = rand_additive


"""
    rand_conflicts_ws98(m; k=2:2:div(m, 2), β=Uniform(), rng=default_rng())
    rand_conflicts_ws98(V::Profile; ...)

Generate a random `Conflicts` contraint, whose underlying graph is constructed
according to the [Watts–Strogatz model](https://doi.org/10.1038/30918). The
keyword arguments `k` and `β` specify the possible values for the corresponding
parameters \$k\$ and \$\\beta\$, which are generated using `rand`. The defaults
are taken from [Hummel and Hetland](https://arxiv.org/abs/2104.06280). Note that
the parameter \$k\$ should be an even number, which Watts and Strogatz assume to
be much smaller than \$m\$.
"""
function rand_conflicts_ws98(m; k=2:2:div(m, 2), β=Uniform(), rng=default_rng())
    G = watts_strogatz(m, rand(rng, k), rand(rng, β), rng=rng)
    return Conflicts(G)
end

rand_conflicts_ws98(V::Profile; kwds...) = rand_conflicts_ws98(ni(V); kwds...)


"""
    rand_conflicts_er59(m, p=Uniform(), rng=default_rng())
    rand_conflicts_er59(V::Profile; ...)
    rand_conflicts(m; ...)
    rand_conflicts(m::Profile; ...)

Generate a random `Conflicts` contraint, whose underlying graph is constructed
according to the [Erdős–Rényi
model](https://doi.org/10.5486%2FPMD.1959.6.3-4.12). The keyword argument `p`
specifies the possible values for the corresponding parameter \$p\$, which is
generated using `rand`.

`rand_conflicts` is an alias for `rand_conflicts_er59`.
"""
function rand_conflicts_er59(m; p=Uniform(), rng=default_rng())
    G = erdos_renyi(m, rand(rng, p), rng=rng)
    return Conflicts(G)
end

rand_conflicts_er59(V::Profile; kwds...) = rand_conflicts_er59(ni(V); kwds...)
const rand_conflicts = rand_conflicts_er59


"""
    rand_conflicts_ba02(m; k=1:m, rng=default_rng())
    rand_conflicts_ba02(V::Profile; ...)

Generate a random `Conflicts` contraint, whose underlying graph is constructed
according to the [Barabási–Albert
model](https://arxiv.org/abs/cond-mat/0106096). The keyword argument `k`
specifies the possible values for the corresponding parameter \$k\$, which is
generated using `rand`.
"""
function rand_conflicts_ba02(m; k=1:m, rng=default_rng())
    G = barabasi_albert(m, rand(rng, k), rng=rng)
    return Conflicts(G)
end

rand_conflicts_ba02(V::Profile; kwds...) = rand_conflicts_ba02(ni(V); kwds...)

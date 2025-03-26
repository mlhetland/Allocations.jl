import JuMP
using Allocations
using Graphs
using Random: seed!, randperm, Xoshiro
using Supposition
using Test

# For utilities tests:
using Allocations: bipartite_matching, lex_optimize!

# For Counts test:
using Allocations: Category

# For symmetry-related tests:
using Allocations: Symmetry, Symmetric, Asymmetric, SymmetrizedConstraint

# For matroid tests:
using Allocations: matroid_partition_knuth73, find_shortest_path,
                   exchange_graph, knuth_matroid, knuth_matroid_erect,
                   matroid_c1, matroid_c2, matroid_c3

# Shorthand for matroid tests:
const btos = Allocations.bits_to_set


function gen_matroid(m, r)
    rand_matroid_knu75(m, r=r:r)
end

itemgen = Data.Integers(1, 10)
rankgen = Data.Integers(1, 10)
matroidgen = @composed gen_matroid(itemgen, rankgen)


function runtests(; slow_tests = true)

    seed!(4252170447285279131)

@testset "Types" begin

    @testset "Allocations" begin

        @testset "Basics" begin

            n, m = 3, 5

            A = Allocation(n, m)

            @test string(A) == "[{}, {}, {}]"

            @test na(A) == length(agents(A)) == n
            @test ni(A) == length(items(A)) == m

            @test isempty(bundle(A, 2))

            give!(A, 2, 4)
            give!(A, 2, 3)
            give!(A, 1, 2)

            @test string(A) == "[{2}, {3, 4}, {}]"

            give!(A, 3, [1, 5])

            @test string(A) == "[{2}, {3, 4}, {1, 5}]"

            @test 4 in bundle(A, 2)

            deny!(A, 2, 4)

            @test string(A) == "[{2}, {3}, {1, 5}]"
            @test summary(A) == "Allocation with 3 agents and 5 items, " *
                                "1 unallocated"

            A₂ = Allocation(A)
            @test A == A₂
            @test A !== A₂

            A₃ = copy(A)
            @test A == A₃
            @test A !== A₃

            deny!(A₂, 2, 3)
            @test A != A₂
            @test A₂ != A₃

            deny!(A₃, 2, 3)
            @test A₂ == A₃

            deny!(A, 2, 3)
            @test A == A₂

            @test isempty(bundle(A, 2))

            @test length(A) == length(collect(A)) == 3

            @test !isempty(A)
            @test isempty(Allocation())

            for g = 1:m
                give!(A, 1, g)
            end

            @test summary(A) == "Allocation with 3 agents and 5 items"

            A = Allocation(2, 4)
            give!(A, 2, 1)
            fill_even!(A)
            @test owner(A, 2) == 1
            @test length(bundle(A, 1)) == length(bundle(A, 2)) == 2

        end

    end

    @testset "Profiles" begin

        let

            X = [1 2 3; 3 2 1]

            V = Additive(X)
            V′ = Profile(X)

            @test V == V′

            i, g, h = 2, 3, 2

            @test value(V, i, g) == 1

            A = Allocation(5, 10)

            @test value(V, i, A) == 0

            give!(A, i, g)
            give!(A, i, h)

            @test value(V, i, A) == 3

        end

        let

            n, m = 5, 10

            V = Additive(n, m)

            i, g = 2, 3

            @test length(agents(V)) == na(V) == n
            @test length(items(V)) == ni(V) == m
            @test value(V, i, g) == 0

            value!(V, i, g, 4)

            @test value(V, i, g) == 4

        end

    end

    @testset "Counts" begin

        C = Counts(
            [1, 2, 3] => 2,
            [4, 5, 6] => 1
        )

        @test C isa Counts
        for c in C
            @test c isa Category
        end

        @test C[2].threshold == 1

    end

    @testset "Conflicts" begin

        C = Conflicts(star_graph(10))

        @test C isa Conflicts

        @test C.graph isa AbstractGraph

    end

    @testset "Inclusion/exclusion" begin

        A = Allocation()

        C₁ = Forbidden(A)
        C₂ = Permitted(A)
        C₃ = Required(A)

        @test C₁ isa Forbidden
        @test C₁.alloc === A
        @test C₂ isa Permitted
        @test C₂.alloc === A
        @test C₃ isa Required
        @test C₃.alloc === A

    end

    @testset "Symmetry" begin

        for C in [Counts, Conflicts]
            @test Symmetry(C) == Symmetric()
        end

        for C in [Forbidden, Permitted, Required]
            @test Symmetry(C) == Asymmetric()
        end

        @test Symmetry(Counts([1] => 2))         == Symmetric()
        @test Symmetry(Conflicts(path_graph(1))) == Symmetric()
        @test Symmetry(Forbidden(Allocation()))  == Asymmetric()
        @test Symmetry(Required(Allocation()))   == Asymmetric()
        @test Symmetry(Permitted(Allocation()))  == Asymmetric()

    end

end

@testset "Data sets" begin

    @test rand_profile === rand_additive

    for _ = 1:10
        V = rand_additive()
        @test V isa Additive
        n, m = na(V), ni(V)
        @test 2 ≤ n ≤ 10
        @test 2n ≤ m ≤ 4n
        @test all(0 .≤ V.values .< 1)
    end

    let V = rand_additive(n=5:5, m=n->10:10)
        @test size(V.values) == (5, 10)
    end

    V = rand_additive()
    m = ni(V)

    for f in [rand_conflicts_ws98, rand_conflicts_er59, rand_conflicts_ba02],
        C in [f(m), f(V)]

        @test C isa Conflicts
        @test C.graph isa AbstractGraph
        @test nv(C.graph) == m

    end

    rng = Xoshiro(1234)
    V₁ = rand_additive(rng=rng)
    rng = Xoshiro(1234)
    V₂ = rand_additive(rng=rng)
    @test V₁.values == V₂.values

end

@testset "Utilities" begin

    @testset "Matching" begin

        X = [0 1 1 0 0 0
             0 0 0 0 0 0
             1 0 0 1 0 0
             0 0 1 0 0 0
             0 0 1 1 0 0
             0 0 0 0 0 1]

        # M = bipartite_matching(Bool.(X))
        M = falses(size(X))
        for (i, g) in bipartite_matching(X)
            M[i, g] = true
        end

        @test all(sum(M, dims=1) .<= 1)
        @test all(sum(M, dims=2) .<= 1)
        @test sum(M) == 5

    end

    @testset "Lexicographic optimization" begin

        model = JuMP.Model(Allocations.conf.MIP_SOLVER)

        JuMP.@variable(model, x >= 0)
        JuMP.@variable(model, y >= 0)

        objectives = [(JuMP.MOI.MIN_SENSE, -x), (JuMP.MOI.MAX_SENSE, y)]

        JuMP.@constraint(model, x + y <= 10)
        JuMP.@constraint(model, x <= 3)

        constraints = lex_optimize!(model, objectives, ϵ=0)

        @test JuMP.termination_status(model) in Allocations.conf.MIP_SUCCESS

        @test JuMP.value(x) ≈ 3
        @test JuMP.value(y) ≈ 7

        # Cleanup
        for con in constraints
            JuMP.delete(model, con)
        end

    end

end

@testset "Basic checks" begin

    A = Allocation(2, 2)

    give!(A, 1, 1)

    @test !check_partition(A)
    @test !check_complete(A)

    give!(A, 2, 2)

    @test check_partition(A)
    @test check_complete(A)

    give!(A, 1, 2)

    @test !check_partition(A)
    @test check_complete(A)

end

@testset "EF checks" begin

    n, m = 2, 3

    V = Additive(ones(n, m))

    A = Allocation(n, m)

    @test check_ef(V, A)

    give!(A, 1, 1)

    @test !check_ef(V, A)
    @test check_ef1(V, A)

    give!(A, 1, 2)
    give!(A, 2, 3)

    @test check_ef1(V, A)
    @test check_efx(V, A)

    value!(V, 2, 1, 2)

    @test check_ef1(V, A)
    @test !check_efx(V, A)

end

@testset "Measures" begin

    V = Additive([1 3; 3 1])
    A = Allocation(2, 2)
    give!(A, 1, 2)
    give!(A, 2, 1)

    @test nash_welfare(V, A) ≈ 9
    @test utility(V, A) ≈ 6
    @test prop_alpha(V, A) ≈ 3 / (4 / 2)

end

@testset "MIPs" begin

    V₀ = Additive(rand(1:10, 3, 15))

    @testset "MNW" begin

        V = V₀

        let res = alloc_mnw(V)

            @test res.alloc isa Allocation
            @test check_ef1(V, res.alloc)
            @test res.mnw > 0

        end

        let res = alloc_mnw([1 2 3; 4 3 1])

            @test string(res.alloc) == "[{3}, {1, 2}]"
            @test res.mnw ≈ 3 * (4 + 3)

            @test res.model isa JuMP.Model

        end

    end

    @testset "MNW with constraints" for C in [
            Conflicts(path_graph(ni(V₀))),
            Counts(
                [1, 2, 3, 4]     => 3,
                [5, 6, 7]        => 2,
                [8, 9, 10]       => 2,
                [11, 12, 13, 14] => 3,
                [15]             => 1
            )
        ]

        V = V₀

        res = alloc_mnw(V)
        resc = alloc_mnw(V, C)

        @test check(V, resc.alloc, C)

        @test resc.alloc isa Allocation
        @test resc.mnw > 0

        # Adding constraint can't improve objective.
        @test resc.mnw <= res.mnw

        @test res.model isa JuMP.Model
        @test resc.model isa JuMP.Model

    end

    @testset "MIP with inclusion/exclusion" begin

        V = Additive([1 1 2; 1 2 1])

        # Somewhat arbitrarily using MNW

        A₀ = Allocation(V)
        give!(A₀, 1, 1)
        A = alloc_mnw(V, Required(A₀)).alloc
        @test 1 in bundle(A, 1)

        A₀ = Allocation(V)
        give!(A₀, 1, 3)
        A = alloc_mnw(V, Forbidden(A₀)).alloc
        @test !(3 in bundle(A, 1))

        A₀ = Allocation(V)
        give!(A₀, 1, [1, 2])
        give!(A₀, 2, [1, 2, 3])
        A = alloc_mnw(V, Permitted(A₀)).alloc
        @test !(3 in bundle(A, 1))

    end

    @testset "MIP with multiple constraints" begin

        V = Profile([1 2 3; 3 1 1])
        F = Forbidden(Allocation(V, 1 => 3))
        R = Required(Allocation(V, 2 => 2))
        C = Constraints(F, R)

        @test string(alloc_mnw(V).alloc) == "[{2, 3}, {1}]"
        @test string(alloc_mnw(V, C).alloc) == "[{1}, {2, 3}]"

    end

    @testset "EF1 with conflicts" begin

        nᵢ = 6
        n, m = nᵢ, 2nᵢ

        V = Additive(rand(n, m))

        # Random graph whose components have at most n items, to ensure EF1 is
        # possible.
        G = blockdiag(SimpleGraph(nᵢ, nᵢ + 1), SimpleGraph(nᵢ, nᵢ + 1))
        C = Conflicts(G)

        res = alloc_ef1(V, C)

        @test res.model isa JuMP.Model

        @test check_ef1(V, res.alloc)

        # Check that it's usable without constraints, even though we're not
        # supplying a single-argument MIP implementation:
        @test check_ef1(V, alloc_ef1(V, nothing).alloc)

    end

    @testset "MNW+EF1 with conflicts" begin

        # Specific example (V and G) where MNW does not lead to EF1:

        V = Additive([10  5  8  1  2  9; 6  1  9  6  8  9])

        G = SimpleGraph(ni(V))
        add_edge!(G, 3, 4)
        add_edge!(G, 5, 6)

        C = Conflicts(G)

        res = alloc_mnw(V, C)

        @test !check_ef1(V, res.alloc)

        res1 = alloc_mnw_ef1(V, C)

        @test res.model isa JuMP.Model
        @test res1.alloc isa Allocation
        @test check_ef1(V, res1.alloc)

        # MNW did not yield EF1, so enforcing EF1 should reduce MNW:
        @test res1.mnw < res.mnw

    end

    @testset "EFX" begin

        V = V₀

        alloc_efx([1 1 1; 1 1 1]) # Bug fix: EFX but not EF

        res = alloc_efx(V)

        A = res.alloc

        @test res.model isa JuMP.Model
        @test A isa Allocation

        @test check_partition(A)
        @test check_complete(A)

        @test check_efx(V, A)

    end

    @testset "Maximin" begin

        V = V₀

        res = alloc_mm(V)

        A = res.alloc
        N = agents(V)

        @test res.model isa JuMP.Model
        @test A isa Allocation
        @test res.mm == minimum(value(V, i, A) for i in N)

    end

    @testset "Lexicographic" begin

        # Sanity check of lexicographic optimization, using internal functions.

        V = Additive([1 0; 1 1])

        ctx = Allocations.init_mip(V, Allocations.conf.MIP_SOLVER)

        A = ctx.alloc_var

        ctx.objectives = [
            (JuMP.MOI.MAX_SENSE, value(V, 1, A))
            (JuMP.MOI.MAX_SENSE, value(V, 2, A))
        ]

        ctx = Allocations.solve_mip(ctx)

        @test string(ctx.alloc) == "[{1}, {2}]"

    end

    @testset "Leximin" begin

        for (case, expected, values) in [

                ([2 1; 1 2]         , "[{1}, {2}]"       , [2, 2])
                ([1 3; 1 2]         , "[{2}, {1}]"       , [1, 3])
                ([4 3; 3 2]         , "[{2}, {1}]"       , [3, 3])
                ([2 2 2 2; 1 1 1 1] , nothing            , [2, 4])
                ([3 3 3 3; 1 1 1 1] , nothing            , [3, 3])
                ([2 1 2 1; 1 2 1 2] , "[{1, 3}, {2, 4}]" , [4, 4])
                ([3 4 2 1; 5 5 1 1] , "[{2, 3}, {1, 4}]" , [6, 6])

                # From Feige, Sapir & Tauber, "A Tight Negative Example for MMS
                # Fair Allocations", WINE'22:
                ([1 16 23 26 4 10 12 19 9
                  1 16 22 26 4  9 13 20 9
                  1 15 23 25 4 10 13 20 9] , nothing     , [39, 40, 43])

            ]

            V = Profile(case)
            N = agents(V)

            res = alloc_lmm(V)
            A = res.alloc

            @test res.model isa JuMP.Model
            @test A isa Allocation

            isnothing(expected) || @test string(A) == expected

            @test sort([value(V, i, A) for i in N]) == values

        end

    end

    slow_tests &&
    @testset "MMS" begin

        let V = Additive([0; 0;;]), res = alloc_mms(V)

            @test res.alloc isa Allocation
            @test mms_alpha(V, res.alloc, res.mmss) == res.alpha == Inf

        end

        let V = Additive([0 0; 1 1]), res = alloc_mms(V)

            @test res.alloc isa Allocation
            @test mms_alpha(V, res.alloc, res.mmss) == res.alpha == 2.0

        end

        V = V₀

        let res = alloc_mms(V)

            @test res.alloc isa Allocation
            @test mms_alpha(V, res.alloc, res.mmss) ≈ res.alpha

        end

        let res = alloc_mms([3 1 2; 4 4 5])

            @test res.model isa JuMP.Model
            @test length(res.mms_models) == 2
            @test all(isa.(res.mms_models, JuMP.Model))

            @test res.alpha ≈ 1.0
            @test res.mmss ≈ [3.0, 5.0]

        end

        let res = alloc_mms([2 1; 1 2])

            @test res.alpha ≈ 2.0

        end

        let res = alloc_mms([2 1; 1 2], cutoff=true)

            @test res.alpha ≈ 1.0

        end

        let

            V = Additive(ones(2, 3))
            C = Conflicts(cycle_graph(3))

            @test_throws "INFEASIBLE" alloc_mms(V, C)

            # No exception:
            @test (alloc_mms(V, C, min_owners=0); true)

            @test_throws "INFEASIBLE" alloc_mms(V, C, min_owners=0,
                                                mms_kwds=[:min_owners=>1])

        end

        let

            V = Additive(ones(2, 1))

            C = Forbidden(Allocation(2, 1, 1 => 1))

            @test_throws "INFEASIBLE" alloc_mms(V, C)

            # No exception:
            @test (alloc_mms(V, C, mms_kwds=(min_owners=0,)); true)

        end

    end

    @testset "GGI" begin

        V = V₀

        let res = alloc_ggi(V)

            @test res.alloc isa Allocation
            @test res.model isa JuMP.Model

        end

        let res = alloc_ggi([1 1 3; 1 1 2])

            @test string(res.alloc) == "[{3}, {1, 2}]"

        end

    end

    @testset "Random with optional constraints" begin

        V = V₀

        let res = alloc_rand_mip(V)

            @test res.alloc isa Allocation
            @test res.model isa JuMP.Model
            @test na(res.alloc) == na(V)
            @test ni(res.alloc) == ni(V)

        end

        C = rand_conflicts_er59(V)

        let res = alloc_rand_mip(V, C, min_owners=0)

            @test res.alloc isa Allocation
            @test res.model isa JuMP.Model
            @test na(res.alloc) == na(V)
            @test ni(res.alloc) == ni(V)

        end

        @test_throws DomainError alloc_rand_mip(Profile(zeros(2, 63)))

    end

    slow_tests &&
    @testset "limits, $func" for func in [
        alloc_ef1, alloc_efx, alloc_mnw, alloc_mnw_ef1, alloc_mm, alloc_ggi,
        alloc_mms, alloc_rand_mip]
        # Could add `cutoff=true` for `alloc_mms`

        V = V₀

        # `alloc_ef1` and `alloc_mnw_ef1` require a constraint, so we'll supply
        # `nothing` to all:
        C = nothing

        N, M = agents(V), items(V)
        n, m = na(V), ni(V)

        for (min_bundle, max_bundle, min_owners, max_owners) in [
                (nothing, 10, 2, nothing)
                (3, 3, nothing, nothing)
                (rand(1:2, n), rand(2:3, n), nothing, rand(1:3, m))
            ]

            lob = something(min_bundle, 0)
            hib = something(max_bundle, m)
            loo = something(min_owners, 0)
            hio = something(max_owners, n)

            A = func(V, C,
                  min_bundle=min_bundle,
                  max_bundle=max_bundle,
                  min_owners=min_owners,
                  max_owners=max_owners).alloc

            broadcast(N, lob, hib) do i, lo, hi
                @test lo <= length(bundle(A, i)) <= hi
            end

            broadcast(M, loo, hio) do i, lo, hi
                @test lo <= length(owners(A, i)) <= hi
            end

        end
    end

end

@testset "Algorithms" begin

    V₀ = Additive(rand(1:10, 4, 12))

    @testset "Random" begin

        V = V₀

        res = alloc_rand(V)

        A = res.alloc

        @test A isa Allocation
        @test check_partition(A)

    end

    @testset "Random with conflicts" begin

        n, m = 10, 50
        nv, ne = m, 30

        V = Additive(zeros(n, m)) # Irrelevant

        # Random graph
        # (Could fail the Δ < n check, if we're unlucky)
        C = Conflicts(SimpleGraph(nv, ne))

        res = alloc_rand(V, C)
        A = res.alloc

        @test A isa Allocation
        # Implied by check_partition, but useful checkpoint:
        @test check_complete(A)
        @test check_partition(A)
        @test check(V, A, C)

    end

    slow_tests &&
    @testset "BKV18(1)" begin

        α = 1.061 # Approximation ratio, for geomean version of NW

        for _ = 1:10

            n = rand(2:4)
            m = 3n
            # Identical valuations
            V = Profile(repeat(rand(1:5, 1, m), n))

            res = alloc_bkv18_1(V)
            A = res.alloc

            @test A isa Allocation
            @test check_partition(A)
            @test α*nash_welfare(V, A)^(1/n) ≥ alloc_mnw(V).mnw^(1/n)
            @test check_efx(V, A)

        end

    end

    @testset "BKV18(2)" begin

        V = Profile([1 1 0 1 0
                     0 1 0 0 1
                     1 1 1 0 1])

        res = alloc_bkv18_2(V)
        A = res.alloc

        @test A isa Allocation
        @test check_partition(A)
        @test res.mnw == alloc_mnw(V).mnw

        V = Profile([1 0; 1 0; 0 1])

        res = alloc_bkv18_2(V)
        A = res.alloc
        @test res.mnw == nash_welfare(V, A) == 1
        @test nash_welfare(V, A, nonzero=false) == 0

        @test alloc_bkv18_2(Profile([0 0; 0 0])).mnw == 0

        V = Profile(ones(2, 10))
        for _ = 1:10
            @test alloc_bkv18_2(V).mnw == 25
        end

        for _ = 1:10
            V = Profile(rand(Bool, rand(2:5), rand(5:15)))
            @test alloc_bkv18_2(V).mnw == alloc_mnw(V).mnw
        end

        # Regression test
        V = Profile([0 0 1 0 1 1; 0 1 0 1 1 1; 0 1 0 1 0 0])
        for _ = 1:10
            alloc_bkv18_2(V) # Should not throw an exception
        end

        # Even distribution of unvalued items:
        V = Profile([1 1 1 0 0 0 0 0 0
                     0 0 0 0 0 0 0 0 0
                     0 0 0 0 0 0 0 0 0])
        A = alloc_bkv18_2(V).alloc
        for g = 4:9
            @test !owned(A, g)
        end

        res = alloc_bkv18_2(V, complete=true)
        @test res.mnw == 3
        A = res.alloc
        for i in agents(A)
            @test length(bundle(A, i)) == 3
        end

    end

    @testset "GHSS18(4)" begin

        # A modified version of the valuation function from the instance
        # where the best possible allocation has α = 3/4. The MMS with this
        # valuation function and m = 2n items, is 2.
        function v(m, B)
            if length(B) != 2
                return length(B)
            end

            g, g′ = minimum(B), maximum(B)
            if g == 1 && g′ == m || floor(g / 2) == floor(g′ / 2)
                return 2
            end

            return 3/2
        end

        # A function that creates a new set, with all items being incremented by
        # 1 (if m ∈ B, then 1 is in the new set). Setting the valuation of an
        # agent to `v(m, inc(m, B))` gives the agent an MMS of 2.
        inc(m, B) = Set(mod(g, 1:m) for g in B)

        n = 5
        m = 2n
        Vf = vcat([B -> v(m, B) for _ in 1:(n - 1)], [B -> v(m, inc(m, B))])
        V = Submodular(Vf, m)

        # Test when the MMS of each agent is known
        res = alloc_ghss18_4(V, repeat([2], n))

        @test !res.fail

        # The value each agent receives should be at least 1/3 * μᵢ = 2/3
        A = res.alloc
        for i in agents(V)
            @test value(V, i, bundle(A, i)) ≥ 2/3
        end

        # Test when the MMS of each agent is unknown
        res = alloc_ghss18_4b(V)

        # The value each agent receives should be at least 1/3 * μᵢ = 2/3
        for i in agents(V)
            @test value(V, i, bundle(A, i)) ≥ 2/3
        end
    end

    slow_tests &&
    @testset "MMS approximation with card. constr." begin

        # A default test set for all algorithms
        V₁ = V₀
        C₁ = Counts(
            [1, 3, 7, 9]          => 1,
            [4, 6, 8, 10, 11, 12] => 3,
            [2, 5]                => 2,
        )
        MMSs₁ = [mms(V₁, i, C₁).mms for i in agents(V₁)]

        # Check if alg finds an allocation that is complete, does not break the
        # cardinality constraints and gives each agent a bundle valued at no
        # less than `α * MMSs[i]`.
        function test_alg(alg, α, V, C, MMSs)

            A = alg(V, C).alloc

            @test A isa Allocation
            @test check_partition(A)

            # The allocation must not break the cardinality constraints
            for i in agents(V), c in C
                @test sum(owner(A, g) == i for g in c) ≤ c.threshold
            end

            for i in agents(V)
                @test value(V, i, bundle(A, i)) ≥ α * MMSs[i]
            end

        end

        @testset "1/3-MMS - BB18(3)" begin

            test_alg(alloc_bb18_3, 1/3, V₁, C₁, MMSs₁)

        end

        @testset "1/2-MMS - HH22(1)" begin

            test_alg(alloc_hh22_1, 1/2, V₁, C₁, MMSs₁)


            # The second half of the bag filling algorithm, where the floor_n(C)
            # highest-valued items are not worth 1/2 and thus ceil_n(C) must be
            # used instead, does not often get run when a random instance is
            # created. Thus, this tests the workings of that part
            let

                C = Counts(
                    [1, 5, 6]           => 3,
                    [2]                 => 1,
                    [3]                 => 1,
                    [4]                 => 3,
                    [7, 8, 9, 10]       => 5,
                )

                V = Additive([
                    0.2 0.4 0.4 0.4 0.1 0.1 0.1 0.1 0.1 0.1;
                    0.2 0.4 0.4 0.4 0.1 0.1 0.1 0.1 0.1 0.1
                ])

                test_alg(alloc_hh22_1, 1/2, V, C, [1, 1])

            end

        end
    end

    slow_tests &&
    @testset "MMS approximation" begin

        @testset "2/3-MMS - GMT18" begin

            function checkvalidallocation(A, V)
                # Test that we have received an allocation
                @test A isa Allocation

                # Test that all items are allocated properly
                for g in items(V)
                    @test owner(A, g) isa Int
                end

                # Test that each agent receives at least (2/3)-MMS
                for i in agents(V)
                    @test value(V, i, bundle(A, i)) ≥ 2/3 * mms(V, i).mms
                end
            end

            @testset "Base test" begin
                V = V₀

                A = alloc_gmt18(V)

                checkvalidallocation(A, V)
            end

            @testset "All bundles contain an item valued at ≥ 2/3" begin
                V = Additive([
                    0.45 0.27 0.10;
                    0.49 0.49 0.49
                ])

                A = alloc_gmt18(V)

                checkvalidallocation(A, V)

                # Test that each agent has an item valued at (2/3)-MMS or higher
                for i in agents(V)
                    @test any(value(V, i, g) ≥ 2/3 * mms(V, i).mms
                                for g in bundle(A, i))
                end
            end

            # Checks issue #1:
            # https://github.com/mlhetland/Allocations.jl/issues/1
            @testset "No item valued at ≥ 1/3" begin
                V = Additive([
                    0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25
                    0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25
                ])

                A = alloc_gmt18(V)

                checkvalidallocation(A, V)
            end
        end

    end

end

@testset "Matroids" begin
    @testset "GraphicMatroid properties" begin
        G = smallgraph(:karate)
        M = GraphicMatroid(G)

        @test rank(M) == 33
        @test rank(M, []) == 0
        @test is_indep(M, [1,2,3])
        @test is_indep(M, [1,2,16])
        @test is_indep(M, [1,2,17]) == false
        @test is_closed(M, 1:50) == false
        @test is_closed(M, 1:78) == true

        # situation encountered during manual testing:
        g = SimpleGraph{Int64}(64, [[9, 10, 11, 12, 13, 14, 15, 16], [9, 10, 11, 12, 13, 14], [9, 10, 11, 12, 13, 14], [9, 10, 11, 16], [9, 10, 11, 12, 13, 14, 15, 16], [9], [9, 10, 12, 14, 15, 16], [9, 10, 11], [1, 2, 3, 4, 5, 6, 7, 8, 10, 11, 12, 13, 14, 15], [1, 2, 3, 4, 5, 7, 8, 9, 11, 12, 13, 14, 15], [1, 2, 3, 4, 5, 8, 9, 10, 12, 13, 14, 15], [1, 2, 3, 5, 7, 9, 10, 11, 13, 15, 16], [1, 2, 3, 5, 9, 10, 11, 12, 16], [1, 2, 3, 5, 7, 9, 10, 11, 15, 16], [1, 5, 7, 9, 10, 11, 12, 14, 16], [1, 4, 5, 7, 12, 13, 14, 15]])
        m = GraphicMatroid(g)
        A = Set([64, 61, 55, 29, 52, 12, 37, 19, 4, 6, 13, 45])

        for e in [33,41,21]
            @test is_indep(m, A ∪ e)
        end
    end

    @testset "UniformMatroid properties" begin
        U = UniformMatroid(10, 6)
        F = FreeMatroid(10)
        Z = ZeroMatroid(10)

        @test is_indep(U, Set())
        @test is_indep(F, Set())
        @test is_indep(Z, Set())
        @test rank(U) == 6
        @test rank(F) == 10
        @test rank(Z) == 0

        for i in 1:5
            S = randperm(10)[1:i]
            @test is_indep(U, S) || S
            @test is_indep(F, S) || S
            @test is_indep(Z, S) == false || S

            @test is_circuit(U, S) == false || S
            @test is_circuit(F, S) == false || S
            @test is_circuit(Z, S) == (i == 1) || S

            @test rank(U, S) == i || S
            @test rank(F, S) == i || S
            @test rank(Z, S) == 0 || S

            @test closure(U, S) == S || S
            @test closure(F, S) == S || S
            @test closure(Z, S) == Set(1:10) || S
        end

        S = randperm(10)[1:6]
        @test is_indep(U, S) == true || S
        @test is_circuit(U, S) == false|| S
        @test rank(U, S) == 6 || S
        @test closure(U, S) == Set(1:10) || S

        S = randperm(10)[1:7]
        @test is_indep(U, S) == false || S
        @test is_circuit(U, S) || S
        @test rank(U, S) == 6 || S
        @test closure(U, S) == Set(1:10) || S


        for i in 8:10
            S = randperm(10)[1:i]
            @test is_indep(U, S) == false || S
            @test is_circuit(U, S) == false || S
            @test rank(U, S) == 6 || S
            @test closure(U, S) == Set(1:10) || S
        end
    end

    @testset "Exchange graphs and transfer paths" begin
        # Every agent likes every item.
        matroids = [FreeMatroid(5) for _ in 1:5]

        # Every agent has 1 item.
        A = Allocation(5,5)
        for i in 1:5
            give!(A, i, i)
        end

        # Every item should be exchangeable with every other item.
        @test exchange_graph(matroids, A) == complete_digraph(5)

        # Zachary's Karate Club.
        k = smallgraph(:karate)
        @test find_shortest_path(k, [1], [16]) |> length == 4
        @test find_shortest_path(k, [1,2,3], [1,18,4]) == [1]

        # A graph with no edges has no transfer paths.
        g = SimpleGraph(10, 0)
        @test find_shortest_path(g, [1,2,3], [4,5,6]) === nothing
        @test find_shortest_path(g, [1,2,3,4], [4,5]) == [4]


        # Three small matroids that require some transfers on the exchange
        # graph to get a partition of the ground set.

        # The ground set:       E = [1       2       3       4       5]
        g1 = SimpleGraph{Int64}(5, [[2, 3], [1, 3], [1, 2], [4],    [5]])
        g2 = SimpleGraph{Int64}(5, [[1],    [3, 4], [2, 4], [2, 3], [5]])
        g3 = SimpleGraph{Int64}(5, [[1],    [2],    [4, 5], [3, 5], [3, 4]])

        ms = [GraphicMatroid(g1), GraphicMatroid(g2), GraphicMatroid(g3)]

        (partition, junk) = matroid_partition_knuth73(ms)
        for (i, set) in enumerate(partition)
            @test is_indep(ms[i], set) || "set $set not indep in ms[$i]"
        end
    end

    @testset "ClosedSetsMatroid properties" begin
        # The example from Knuth (1975) section 3.
        n = 10                                    # 134   159    256   358    379    238
        enlargements = [nothing, [btos(x) for x in [0x1a, 0x222, 0x64, 0x128, 0x288, 0x10c]]]
        M = knuth_matroid(n, enlargements)

        B = Set([btos(x) for x in [0x00a9, 0x0131, 0x0207, 0x00c3, 0x008d, 0x0123, 0x00c5, 0x0087, 0x02a1, 0x01c1, 0x0185, 0x0291, 0x0099, 0x0113, 0x0151, 0x001d, 0x0213, 0x02c1, 0x010b, 0x008b, 0x020d, 0x0189, 0x0039, 0x00e1, 0x0285, 0x01a1, 0x0311, 0x0191, 0x0245, 0x0381, 0x0107, 0x00a3, 0x00d1, 0x0143, 0x0231, 0x00b1, 0x004b, 0x0243, 0x0017, 0x0035, 0x0053, 0x0115, 0x0261, 0x002b, 0x020b, 0x00c9, 0x0063, 0x0055, 0x0059, 0x0119, 0x0249, 0x0183, 0x0229, 0x0251, 0x0027, 0x0033, 0x0309, 0x0225, 0x0219, 0x0321, 0x0095, 0x0047, 0x0283, 0x0071, 0x000f, 0x0305, 0x0341, 0x0215, 0x0093, 0x00a5, 0x0303]])

        C = Set([btos(x) for x in [0x001a, 0x0222, 0x002c, 0x004c, 0x010c, 0x0064, 0x0124, 0x0144, 0x0068, 0x0128, 0x0148, 0x0288, 0x0160, 0x008e, 0x020e, 0x0036, 0x0056, 0x0096, 0x0116, 0x0216, 0x00a6, 0x00c6, 0x0246, 0x0186, 0x0286, 0x0306, 0x00aa, 0x00ca, 0x024a, 0x018a, 0x030a, 0x0072, 0x00b2, 0x0132, 0x00d2, 0x0152, 0x0252, 0x0192, 0x0292, 0x0312, 0x00e2, 0x01a2, 0x01c2, 0x02c2, 0x0342, 0x0382, 0x009c, 0x021c, 0x00b4, 0x0234, 0x00d4, 0x0254, 0x0194, 0x0294, 0x0314, 0x02a4, 0x02c4, 0x0384, 0x00b8, 0x0238, 0x00d8, 0x0258, 0x0198, 0x0318, 0x00f0, 0x0270, 0x01b0, 0x02b0, 0x0330, 0x01d0, 0x02d0, 0x0350, 0x0390, 0x02e0, 0x03a0, 0x03c0]])

        @test rank(M) == 4
        @test rank(M, btos(0b01101100)) == 2

        for enl in [0x1a, 0x222, 0x64, 0x128, 0x288, 0x10c]
            @test is_indep(M, btos(enl)) == false
        end

        for base in B
            @test is_indep(M, base) == true
        end
        for circ in C
            @test is_indep(M, circ) == false
        end

        for base in B
            @test rank(M, base) == 4
        end

        @test closure(M, btos(0x1a)) == btos(0x1a)
        @test closure(M, btos(0x222)) == btos(0x222)
        @test closure(M, btos(0x288)) == btos(0x288)

        for x in [0x64, 0x128, 0x10c]
            x = btos(x)
            @test (x, closure(M, x)) == (x, btos(0x016c))
        end

        for base in B
            @test is_circuit(M, base) == false
        end
        for circ in C
            @test is_circuit(M, circ) == true
        end

        @test minimal_spanning_subset(M, SmallBitSet(1:n)) in B
        @test minimal_spanning_subsets(M, SmallBitSet(1:n)) == B
        @test bases(M) == B
    end

    @testset "FullMatroid properties" begin
        # The example from Knuth (1975) section 3.
        n = 10                  # 134   159    256   358    379    238
        enlargements = [nothing, [btos(x) for x in [0x1a, 0x222, 0x64, 0x128, 0x288, 0x10c]]]
        M = knuth_matroid_erect(n, enlargements)

        B = Set([btos(x) for x in [0x00a9, 0x0131, 0x0207, 0x00c3, 0x008d, 0x0123, 0x00c5, 0x0087, 0x02a1, 0x01c1, 0x0185, 0x0291, 0x0099, 0x0113, 0x0151, 0x001d, 0x0213, 0x02c1, 0x010b, 0x008b, 0x020d, 0x0189, 0x0039, 0x00e1, 0x0285, 0x01a1, 0x0311, 0x0191, 0x0245, 0x0381, 0x0107, 0x00a3, 0x00d1, 0x0143, 0x0231, 0x00b1, 0x004b, 0x0243, 0x0017, 0x0035, 0x0053, 0x0115, 0x0261, 0x002b, 0x020b, 0x00c9, 0x0063, 0x0055, 0x0059, 0x0119, 0x0249, 0x0183, 0x0229, 0x0251, 0x0027, 0x0033, 0x0309, 0x0225, 0x0219, 0x0321, 0x0095, 0x0047, 0x0283, 0x0071, 0x000f, 0x0305, 0x0341, 0x0215, 0x0093, 0x00a5, 0x0303]])

        C = Set([btos(x) for x in [0x001a, 0x0222, 0x002c, 0x004c, 0x010c, 0x0064, 0x0124, 0x0144, 0x0068, 0x0128, 0x0148, 0x0288, 0x0160, 0x008e, 0x020e, 0x0036, 0x0056, 0x0096, 0x0116, 0x0216, 0x00a6, 0x00c6, 0x0246, 0x0186, 0x0286, 0x0306, 0x00aa, 0x00ca, 0x024a, 0x018a, 0x030a, 0x0072, 0x00b2, 0x0132, 0x00d2, 0x0152, 0x0252, 0x0192, 0x0292, 0x0312, 0x00e2, 0x01a2, 0x01c2, 0x02c2, 0x0342, 0x0382, 0x009c, 0x021c, 0x00b4, 0x0234, 0x00d4, 0x0254, 0x0194, 0x0294, 0x0314, 0x02a4, 0x02c4, 0x0384, 0x00b8, 0x0238, 0x00d8, 0x0258, 0x0198, 0x0318, 0x00f0, 0x0270, 0x01b0, 0x02b0, 0x0330, 0x01d0, 0x02d0, 0x0350, 0x0390, 0x02e0, 0x03a0, 0x03c0]])

        @test rank(M) == 4
        @test rank(M, btos(0b01101100)) == 2

        for enl in [0x1a, 0x222, 0x64, 0x128, 0x288, 0x10c]
            @test is_indep(M, btos(enl)) == false
        end

        for base in B
            @test is_indep(M, base) == true
        end
        for circ in C
            @test is_indep(M, circ) == false
        end

        for base in B
            @test rank(M, base) == 4
        end

        @test closure(M, btos(0x1a)) == btos(0x1a)
        @test closure(M, btos(0x222)) == btos(0x222)
        @test closure(M, btos(0x288)) == btos(0x288)

        for x in [0x64, 0x128, 0x10c]
            x = btos(x)
            @test (x, closure(M, x)) == (x, btos(0x016c))
        end

        for base in B
            @test is_circuit(M, base) == false
        end
        for circ in C
            @test is_circuit(M, circ) == true
        end

        @test minimal_spanning_subset(M, SmallBitSet(1:n)) in B
        @test minimal_spanning_subsets(M, SmallBitSet(1:n)) == B
        @test bases(M) == B
    end

    @testset "knuth_matroid_erect" begin
        # The example from Knuth (1975) section 3.
        n = 10
        enlargements = [nothing, [btos(x) for x in [0x1a, 0x222, 0x64, 0x128, 0x288, 0x10c]]]

        F0 = Set((SmallBitSet(),)) # r=0: The empty set alone.
        F1 = Set([Set(i) for i in 1:n]) # r=1: Singleton subsets of E.
        F2 = Set([Set(x) for x in [[1, 2], [1, 3], [1, 4], [1, 5], [1, 6], [1, 7], [1, 8], [1, 9],
            [1, 10], [2, 3], [2, 4, 5], [2, 6, 10], [2, 7], [2, 8], [2, 9],
            [3, 4, 6, 7, 9], [3, 5], [3, 8], [3, 10], [4, 8, 10], [5, 6],
            [5, 7], [5, 8], [5, 9], [5, 10], [6, 8], [7, 8], [7, 10],
            [8, 9], [9, 10]]])
        F3 = Set([Set(x) for x in
                    [[1, 2, 3], [1, 2, 4, 5], [1, 2, 6, 10], [1, 2, 7], [1, 2, 8],
            [1, 2, 9], [1, 3, 4, 6, 7, 9], [1, 3, 5], [1, 3, 8], [1, 3, 10],
            [1, 4, 8, 10], [1, 5, 6], [1, 5, 7], [1, 5, 8], [1, 5, 9],
            [1, 5, 10], [1, 6, 8], [1, 7, 8], [1, 7, 10], [1, 8, 9],
            [1, 9, 10], [2, 3, 4, 5, 6, 7, 8, 9, 10]]])
        F4 = Set((SmallBitSet(1:n),)) # r=4: The family of only one set - E
        F = [F0, F1, F2, F3, F4]
        B = Set([btos(x) for x in [0x00a9, 0x0131, 0x0207, 0x00c3, 0x008d, 0x0123, 0x00c5, 0x0087, 0x02a1, 0x01c1, 0x0185, 0x0291, 0x0099, 0x0113, 0x0151, 0x001d, 0x0213, 0x02c1, 0x010b, 0x008b, 0x020d, 0x0189, 0x0039, 0x00e1, 0x0285, 0x01a1, 0x0311, 0x0191, 0x0245, 0x0381, 0x0107, 0x00a3, 0x00d1, 0x0143, 0x0231, 0x00b1, 0x004b, 0x0243, 0x0017, 0x0035, 0x0053, 0x0115, 0x0261, 0x002b, 0x020b, 0x00c9, 0x0063, 0x0055, 0x0059, 0x0119, 0x0249, 0x0183, 0x0229, 0x0251, 0x0027, 0x0033, 0x0309, 0x0225, 0x0219, 0x0321, 0x0095, 0x0047, 0x0283, 0x0071, 0x000f, 0x0305, 0x0341, 0x0215, 0x0093, 0x00a5, 0x0303]])
        C = Set([btos(x) for x in [0x001a, 0x0222, 0x002c, 0x004c, 0x010c, 0x0064, 0x0124, 0x0144, 0x0068, 0x0128, 0x0148, 0x0288, 0x0160, 0x008e, 0x020e, 0x0036, 0x0056, 0x0096, 0x0116, 0x0216, 0x00a6, 0x00c6, 0x0246, 0x0186, 0x0286, 0x0306, 0x00aa, 0x00ca, 0x024a, 0x018a, 0x030a, 0x0072, 0x00b2, 0x0132, 0x00d2, 0x0152, 0x0252, 0x0192, 0x0292, 0x0312, 0x00e2, 0x01a2, 0x01c2, 0x02c2, 0x0342, 0x0382, 0x009c, 0x021c, 0x00b4, 0x0234, 0x00d4, 0x0254, 0x0194, 0x0294, 0x0314, 0x02a4, 0x02c4, 0x0384, 0x00b8, 0x0238, 0x00d8, 0x0258, 0x0198, 0x0318, 0x00f0, 0x0270, 0x01b0, 0x02b0, 0x0330, 0x01d0, 0x02d0, 0x0350, 0x0390, 0x02e0, 0x03a0, 0x03c0]])

        M = knuth_matroid_erect(n, enlargements)

        @test M.F[1] == F0
        @test M.F[2] == F1
        @test M.F[3] == F2
        @test M.F[4] == F3
        @test M.F[5] == F4

        @test M.I[5] == B
        @test M.C == C
    end

    @testset "Pi-based KMC example" begin
        # The example from Knuth (1975) section 3.
        n = 10
        enlargements = [nothing, [btos(x) for x in [0x1a, 0x222, 0x64, 0x128, 0x288, 0x10c]]]

        F0 = Set((SmallBitSet(),)) # r=0: The empty set alone.
        F1 = Set([Set(i) for i in 1:n]) # r=1: Singleton subsets of E.
        F2 = Set([Set(x) for x in
                    [[1, 2], [1, 3], [1, 4], [1, 5], [1, 6], [1, 7], [1, 8], [1, 9],
            [1, 10], [2, 3], [2, 4, 5], [2, 6, 10], [2, 7], [2, 8], [2, 9],
            [3, 4, 6, 7, 9], [3, 5], [3, 8], [3, 10], [4, 8, 10], [5, 6],
            [5, 7], [5, 8], [5, 9], [5, 10], [6, 8], [7, 8], [7, 10],
            [8, 9], [9, 10]]])
        F3 = Set([Set(x) for x in
                    [[1, 2, 3], [1, 2, 4, 5], [1, 2, 6, 10], [1, 2, 7], [1, 2, 8],
            [1, 2, 9], [1, 3, 4, 6, 7, 9], [1, 3, 5], [1, 3, 8], [1, 3, 10],
            [1, 4, 8, 10], [1, 5, 6], [1, 5, 7], [1, 5, 8], [1, 5, 9],
            [1, 5, 10], [1, 6, 8], [1, 7, 8], [1, 7, 10], [1, 8, 9],
            [1, 9, 10], [2, 3, 4, 5, 6, 7, 8, 9, 10]]])
        F4 = Set((SmallBitSet(1:n),)) # r=4: The family of only one set - E
        F = [F0, F1, F2, F3, F4]

        M = knuth_matroid(n, enlargements)
        result = M.F

        @test result[1] == F0
        @test result[2] == F1
        @test result[3] == F2
        @test result[4] == F3
        @test result[5] == F4
    end

    @testset "Random matroid properties" begin
        max_ex = 1_000
        @check max_examples=max_ex matroid_c1(matroidgen)
        @check max_examples=max_ex matroid_c2(matroidgen)
        @check max_examples=max_ex matroid_c3(matroidgen)
    end

end

return nothing

end

import JuMP
using Allocations
using Graphs
using Random: seed!
using Test

# For matching test:
using Allocations: bipartite_matching

# For Counts test:
using Allocations: Category

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

            deny!(A, 2, 3)

            @test isempty(bundle(A, 2))

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

        @test graph(C) isa AbstractGraph

    end

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

    slow_tests &&
    @testset "MMS" begin

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
                @test sum(owner(A, g) == i for g in c) ≤ threshold(c)
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
            V = V₀

            A = alloc_gmt18(V)

            @test A isa Allocation
            # Test that all items are allocated properly
            for g in items(V)
                @test owner(A, g) isa Int
            end

            for i in agents(V)
                @test value(V, i, bundle(A, i)) ≥ 2/3 * mms(V, i).mms
            end
        end

    end

end

return nothing

end

import JuMP
using Allocations
using LightGraphs
using Random: seed!
using Test

# For matching test:
using Allocations: bipartite_matching

# For Counts test:
using Allocations: Category

function runtests()

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

        end

    end

    @testset "Valuations" begin

        let

            X = [1 2 3; 3 2 1]

            V = Additive(X)
            V′ = Valuation(X)

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

        M = bipartite_matching(Bool.(X))

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

    V = Additive([1 2; 2 1])
    A = Allocation(2, 2)
    give!(A, 1, 2)
    give!(A, 2, 1)

    @test nash_welfare(V, A) ≈ 4
    @test prop_alpha(V, A) ≈ 2 / (3 / 2)

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

        n, m = 3, 9
        nv, ne = m, 2m

        V = Additive(rand(n, m))

        # Random graph
        C = Conflicts(SimpleGraph(nv, ne))

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

    @testset "MGG" begin

        V = V₀

        let res = alloc_mgg(V)

            @test res.alloc isa Allocation
            @test res.model isa JuMP.Model

        end

        let res = alloc_mgg([1 1 3; 1 1 2])

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

    @testset "1/2-MMS with card. constr." begin

        V = V₀

        function test_cardinality_constraints_half_mms(V, C)

            A = alloc_half_mms(V, C)

            @test A isa Allocation
            # Test that all items are allocated properly
            for g in items(V)
                @test owner(A, g) isa Int
            end

            # The allocation must not break the cardinality constraints
            for i in agents(V), c in C
                @test sum(owner(A, g) == i for g in c) <= threshold(c)
            end

            for i in agents(V)
                @test value(V, i, bundle(A, i)) >= 0.5 * mms(V, i, C).mms
            end

        end

        let

            C = Counts(
                [1, 3, 7, 9]          => 1,
                [4, 6, 8, 10, 11, 12] => 3,
                [2, 5]                => 2,
            )

            test_cardinality_constraints_half_mms(V, C)

        end

        # The second half of the bag filling algorithm, where the floor_n(C)
        # highest-valued items are not worth 1/2 and thus ceil_n(C) must be
        # used instead, does not often get ran when a random instance is
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

            test_cardinality_constraints_half_mms(V, C)

        end

    end

end

return nothing

end

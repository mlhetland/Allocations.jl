using Allocations
using LightGraphs
using Test
using Random: seed!

# For matching test:
using Allocations: bipartite_matching

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

            @test 4 in bundle(A, 2)


            deny!(A, 2, 4)

            @test string(A) == "[{2}, {3}, {}]"
            @test summary(A) == "Allocation with 3 agents and 5 items, " *
                                "3 unallocated"

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

            @test value(V, i, bundle(A, i)) == 0

            give!(A, i, g)
            give!(A, i, h)

            @test value(V, i, bundle(A, i)) == 3

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

@testset "MIPs" begin

    n, m = 3, 15
    X = rand(1:10, n, m)
    V = Additive(X)

    @testset "MNW" begin

        res = alloc_mnw(V)

        @test res.alloc isa Allocation
        @test check_ef1(V, res.alloc)
        @test res.mnw > 0

        res = alloc_mnw([1 2 3; 4 3 1])

        @test string(res.alloc) == "[{3}, {1, 2}]"
        @test res.mnw ≈ 3 * (4 + 3)

    end

    @testset "MNW with constraints" begin

        G = path_graph(m)
        C = Conflicts(G)

        res = alloc_mnw(V)
        resc = alloc_mnw(V, C)

        @test check(V, resc.alloc, C)

        @test resc.alloc isa Allocation
        @test resc.mnw > 0

        # Adding constraint can't improve objective.
        @test resc.mnw <= res.mnw

    end

    @testset "EF1 with conflicts" begin

        n, m = 3, 9
        nv, ne = m, 2m

        V = Additive(rand(n, m))

        # Random graph
        C = Conflicts(SimpleGraph(nv, ne))

        @test check_ef1(V, alloc_ef1(V, C).alloc)

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

        @test res1.alloc isa Allocation
        @test check_ef1(V, res1.alloc)

        # MNW did not yield EF1, so enforcing EF1 should reduce MNW:
        @test res1.mnw < res.mnw

    end

    @testset "Maximin" begin

        res = alloc_mm(V)

        A = res.alloc
        N = agents(V)

        @test A isa Allocation
        @test res.mm == minimum(value(V, i, bundle(A, i)) for i in N)

    end

    @testset "MMS" begin

        res = alloc_mms(V)

        @test res.alloc isa Allocation

        res = alloc_mms([3 1 2; 4 4 5])

        @test res.alpha ≈ 1.0
        @test res.mmss ≈ [3.0, 5.0]

        res = alloc_mms([2 1; 1 2])

        @test res.alpha ≈ 2.0

        res = alloc_mms([2 1; 1 2], cutoff=true)

        @test res.alpha ≈ 1.0

    end

    @testset "MGG" begin

        res = alloc_mgg(V)

        @test res.alloc isa Allocation

        res = alloc_mgg([1 1 3; 1 1 2])

        @test string(res.alloc) == "[{3}, {1, 2}]"

    end

end

return nothing

end

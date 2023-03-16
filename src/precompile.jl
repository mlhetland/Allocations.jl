using SnoopPrecompile

@precompile_setup begin

    __init__() # from conf.jl

    n, m = 2, 2

    @precompile_all_calls begin

        for T in [Int, Float64]
            V = Profile(ones(T, n, m))

            alloc_bkv18_1(V) # Identical
            alloc_bkv18_2(V) # Binary

            # Calls alloc_ghss18_4:
            alloc_ghss18_4b(Submodular([S -> T(length(S)) for i = 1:n], m))

            alloc_gmt18(V)

            alloc_rand(n, m)

            C₁ = Conflicts(SimpleGraph(n, n-1))

            alloc_rand(n, m, C₁)

            C₂ = Counts(collect(1:n) => n)

            alloc_bb18_3(V, C₂)
            A = alloc_half_mms(V, C₂).alloc

            check(V, A, C₂)
            check_complete(A)
            check_partition(A)
            check_ef1(V, A)
            check_efx(V, A)
            check_ef(V, A)
            utility(V, A)

            for C in [nothing, C₁, C₂]
                alloc_ef1(V, C)
                alloc_efx(V, C)
                alloc_ggi(V, C)
                alloc_mms(V, C)
                alloc_mnw(V, C)
                alloc_mnw_ef1(V, C)
            end

        end

    end

end

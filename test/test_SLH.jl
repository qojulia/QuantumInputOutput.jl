using QuantumInputOutput
using SecondQuantizedAlgebra
using QuantumOptics
using QuantumOpticsBase: dagger
using SymbolicUtils
using FunctionWrappers: FunctionWrapper
using Test

@testset "SLH" begin
    hu1 = FockSpace(:u1)
    hc1 = FockSpace(:c1)
    hv1 = FockSpace(:v1)
    h = hu1 ⊗ hc1 ⊗ hv1

    au = Destroy(h, :a_u, 1)
    c = Destroy(h, :c, 2)
    av = Destroy(h, :a_v, 3)

    @variables gu::Real Δ::Real γ::Real
    @variables gv::Complex

    G_u = SLH(1, gu'*au, 0) # input cavity
    G_c = SLH(1, √(γ)*c, Δ*c'c) # system cavity
    G_v = SLH(1, gv'*av, 0) # output cavity

    @test size(scattering(G_c)) == (1, 1)
    @test length(lindblad(G_c)) == 1
    @test iszero(simplify(hamiltonian(G_c) - Δ*c'c))

    SLH(1, [√(γ)*c], Δ*c'c)
    @test isequal(G_c, SLH(1, [√(γ)*c], Δ*c'c))

    @testset "simple_cascade" begin
        G1 = G_u ▷ G_c
        @test size(scattering(G1)) == (1, 1)
        @test iszero(simplify(lindblad(G1)[1] - simplify(gu'*au + √(γ)*c)))
        expected_H = simplify(
            hamiltonian(G_c) - 1im/2*((√(γ)*c)'*(1)*gu'*au - (gu'*au)'*(1)*(√(γ)*c)),
        )
        @test iszero(simplify(hamiltonian(G1) - expected_H))

        G2 = cascade(G_u, G_c, G_v)
        G3 = G_u ▷ G_c ▷ G_v

        @test isequal(G2, G1 ▷ G_v)
        @test isequal(G2, ▷(G1, G_v))
        @test isequal(G2, G3)

        @test iszero(simplify(lindblad(G2)[1] - (gu'*au + √(γ)*c + gv'*av)))
    end

    @testset "simple_concatenate" begin
        G1 = SLH(1, gu'*au, 0)
        G2 = SLH(1, √(γ)*c, Δ*c'c)

        Gc = concatenate(G1, G2)

        @test size(scattering(Gc)) == (2, 2)
        @test length(lindblad(Gc)) == 2
        @test isequal(lindblad(Gc)[1], gu'*au)
        @test isequal(lindblad(Gc)[2], √(γ)*c)
        @test isequal(hamiltonian(Gc), Δ*c'c)

        Gc2 = G1 ⊞ G2
        @test isequal(Gc, Gc2)
    end

    @testset "two-port symbolic cascade" begin
        hu2 = FockSpace(:u2)
        hv2 = FockSpace(:v2)
        h2 = hu1 ⊗ hu2 ⊗ hv1 ⊗ hv2

        au1 = Destroy(h2, :a_u1, 1)
        au2 = Destroy(h2, :a_u2, 2)
        av1 = Destroy(h2, :a_v1, 3)
        av2 = Destroy(h2, :a_v2, 4)

        @variables gu1::Real gu2::Real gv1::Real gv2::Real t::Real r::Real

        G_bs = SLH([r t; t -r], [0, 0], 0)
        G_out = SLH(1, gv1' * av1, 0) ⊞ SLH(1, gv2' * av2, 0)
        G_cas = G_bs ▷ G_out

        @test length(lindblad(G_cas)) == 2
        @test iszero(simplify(lindblad(G_cas)[1] - (gv1' * av1)))
        @test iszero(simplify(lindblad(G_cas)[2] - (gv2' * av2)))
    end

    @testset "show" begin
        @test sprint(show, G_u) == "SLH{1} with 1 port"
        @test sprint(show, MIME("text/plain"), G_u) ==
              "SLH{1} with 1 port\n" *
              "  S = 1×1 identity\n" *
              "  L = 1 symbolic jump operator\n" *
              "  H = 0"

        @test sprint(show, MIME("text/plain"), G_u ⊞ G_c) ==
              "SLH{2} with 2 ports\n" *
              "  S = 2×2 identity\n" *
              "  L = 2 symbolic jump operators\n" *
              "  H = symbolic, 1 term"

        # no equation ever reaches the output, whatever the system size
        M = 6
        @variables s[1:M, 1:M]::Complex
        ops = [gu * au, √(γ) * c, gv * av, gu * au * c' * av, c' * c * av, au * av' * c]
        G_big = SLH([s[i, j] for i = 1:M, j = 1:M], ops, sum(ops[i]' * ops[i] for i = 1:M))
        @test sprint(show, MIME("text/plain"), G_big) ==
              "SLH{6} with 6 ports\n" *
              "  S = 6×6 symbolic\n" *
              "  L = 6 symbolic jump operators\n" *
              "  H = symbolic, 8 terms"

        # passive component: nothing couples to the ports
        @variables r::Real τ::Real
        @test sprint(show, MIME("text/plain"), SLH([r τ; τ -r], [0, 0], 0)) ==
              "SLH{2} with 2 ports\n" *
              "  S = 2×2 symbolic\n" *
              "  L = 2 jump operators (all zero)\n" *
              "  H = 0"

        @test sprint(show, MIME("text/plain"), SLH(1, [0], Δ)) ==
              "SLH{1} with 1 port\n" *
              "  S = 1×1 identity\n" *
              "  L = 1 jump operator (all zero)\n" *
              "  H = symbolic"

        @test sprint(show, MIME("text/plain"), SLH(1, [0], 2)) ==
              "SLH{1} with 1 port\n" *
              "  S = 1×1 identity\n" *
              "  L = 1 jump operator (all zero)\n" *
              "  H = 2"

        L_f = FunctionWrapper{typeof(gu' * au),Tuple{Float64}}(t -> gu' * au)
        @test sprint(show, MIME("text/plain"), SLH(1, L_f, 0 * au)) ==
              "SLH{1} with 1 port\n" *
              "  S = 1×1 identity\n" *
              "  L = 1 time-dependent jump operator\n" *
              "  H = time-dependent"

        bc = FockBasis(4)
        a_op = destroy(bc)
        @test sprint(
            show,
            MIME("text/plain"),
            SLH(1, sparse(a_op), sparse(dagger(a_op) * a_op)),
        ) ==
              "SLH{1} with 1 port\n" *
              "  S = 1×1 identity\n" *
              "  L = 1 numeric jump operator (5×5)\n" *
              "  H = numeric, 5×5"
    end

    @testset "numeric operators" begin
        bc = FockBasis(4)
        a_op = destroy(bc)
        H_s = sparse(0.5 * dagger(a_op) * a_op)
        L_s = sparse(sqrt(1.0) * a_op)
        gu_f(t) = exp(-t^2) * sparse(a_op)
        gv_f(t) = exp(-(t - 2)^2) * sparse(a_op)

        @testset "cascade evaluates time-dependent operators" begin
            G1 = SLH(1, gu_f, H_s)
            G2 = SLH(1, gv_f, H_s)
            G_cas = G1 ▷ G2
            @test lindblad(G_cas)[1](0.5) == gu_f(0.5) + gv_f(0.5)
        end

        @testset "cascade evaluates mixed static/time-dependent operators" begin
            G_cas = SLH(1, L_s, H_s) ▷ SLH(1, gu_f, H_s)
            @test lindblad(G_cas)[1](0.5) == L_s + gu_f(0.5)
        end

        @testset "concatenation evaluates mixed static/time-dependent operators" begin
            G_cat = SLH(1, L_s, H_s) ⊞ SLH(1, gu_f, H_s)
            @test lindblad(G_cat)[1](0.5) == L_s
            @test lindblad(G_cat)[2](0.5) == gu_f(0.5)
        end

        @testset "concatenation preserves static operators" begin
            G_cat = SLH(1, L_s, H_s) ⊞ SLH(1, L_s, H_s)
            @test lindblad(G_cat)[1] == L_s
            @test lindblad(G_cat)[2] == L_s
        end

        @testset "public accessors preserve time-dependent operators" begin
            G_td = SLH(1, gu_f, H_s)
            @test lindblad(G_td)[1](0.5) == gu_f(0.5)
            @test hamiltonian(G_td)(0.5) == H_s
        end

        @testset "public accessors preserve static operators" begin
            G_s = SLH(1, L_s, H_s)
            @test lindblad(G_s)[1] == L_s
            @test hamiltonian(G_s) == H_s
        end

        @testset "SLH with only plain closures errors" begin
            bare_f(t) = t * ones(ComplexF64, 5, 5)
            bare_g(t) = (1 - t) * ones(ComplexF64, 5, 5)
            @test_throws ErrorException SLH([1 0; 0 1], [bare_f, bare_g], bare_g)
        end

    end
end

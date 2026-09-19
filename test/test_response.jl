using QuantumInputOutput
using SecondQuantizedAlgebra
using QuantumOptics
using QuantumOpticsBase: dagger, identityoperator
using LinearAlgebra
using Test

@testset "stationary frequency response" begin
    h = FockSpace(:response_mode)
    a = Destroy(h, :a)
    b = FockBasis(20)

    @variables Δ::Real κ1::Real κ2::Real ϕ::Real Gp::Real

    Hlin = -Δ * a' * a
    Glin = SLH([1 0; 0 1], [√(κ1) * a, √(κ2) * a], Hlin)

    function steady(G, p)
        Hn, Jn = to_numeric(G, b; parameter = p)
        return steadystate.eigenvector(Hn, collect(Jn))
    end

    @testset "output field uses the package SLH sign convention" begin
        p = Dict(Δ => 0.0, κ1 => 0.4, κ2 => 0.6)
        ρ = steady(Glin, p)
        R = frequency_response(Glin, b, ρ; parameter = p)

        L2 = to_numeric(√(κ2) * a, b; parameter = p)
        bout = to_numeric(output_field(Glin, 2), b; parameter = p)
        @test Matrix(bout.data) ≈ Matrix(L2.data)

        bout_coherent = output_field(R, 2; input = [0.0, 2.0])
        expected = L2 + 2.0 * identityoperator(b)
        @test Matrix(bout_coherent.data) ≈ Matrix(expected.data)
    end

    @testset "linear cavity: analytic S parameters and susceptibility" begin
        Δ_ = 0.7
        κ1_ = 0.4
        κ2_ = 0.6
        κ_ = κ1_ + κ2_
        p = Dict(Δ => Δ_, κ1 => κ1_, κ2 => κ2_)
        ρ = steady(Glin, p)
        R = frequency_response(Glin, b, ρ; parameter = p)
        ω = collect(-2.0:0.1:2.0)

        S21 = scattering_parameter(R, ω; in_port = 1, out_port = 2)
        S11 = scattering_parameter(R, ω; in_port = 1, out_port = 1)
        analytic_S21(w) = -sqrt(κ1_ * κ2_) / (κ_ / 2 - im * (Δ_ + w))
        analytic_S11(w) = 1 - κ1_ / (κ_ / 2 - im * (Δ_ + w))

        @test all(isapprox.(S21, analytic_S21.(ω); atol = 1e-7))
        @test all(isapprox.(S11, analytic_S11.(ω); atol = 1e-7))
        @test maximum(abs.(abs2.(S11) .+ abs2.(S21) .- 1)) < 1e-8

        a_num = to_numeric(a, b)
        adag_num = dagger(a_num)
        χ = susceptibility(R, a_num, adag_num, ω)
        analytic_χ(w) = im / (-κ_ / 2 + im * (Δ_ + w))
        @test all(isapprox.(χ, analytic_χ.(ω); atol = 1e-7))

        # Symbolic observables use the parameter/operator context cached in R.
        χ_symbolic = susceptibility(R, a, a', 0.25)
        @test χ_symbolic ≈ analytic_χ(0.25) atol = 1e-7

        # Scalar and vector frequency APIs are deliberately consistent.
        @test scattering_parameter(R, 0.0; in_port = 1, out_port = 2) ≈ analytic_S21(0.0)
        @test susceptibility(R, a, a', 0.0) ≈ analytic_χ(0.0)
        @test length(scattering_parameter(R, ω; in_port = 1, out_port = 2)) == length(ω)

        full = scattering_response(R, ω)
        @test size(full.normal) == (2, 2, length(ω))
        @test size(full.anomalous) == (2, 2, length(ω))
        @test maximum(abs.(full.anomalous)) < 1e-10

        # Prepared and one-shot interfaces are equivalent.
        one_shot =
            scattering_parameter(Glin, b, ρ, ω; parameter = p, in_port = 1, out_port = 2)
        @test one_shot ≈ S21
    end

    @testset "one-port reflection defaults and symbolic scattering matrix" begin
        κ_ = 1.0
        p1 = Dict(Δ => 0.0, κ1 => κ_)
        Gone = SLH(1, √(κ1) * a, -Δ * a' * a)
        ρ1 = steady(Gone, p1)
        R1 = frequency_response(Gone, b, ρ1; parameter = p1)
        @test scattering_parameter(R1, 0.0) ≈ -1 atol = 1e-8

        phase = 0.37
        Gphase = SLH(SecondQuantizedAlgebra.expim(ϕ), √(κ1) * a, -Δ * a' * a)
        @test scattering(Gphase)[1, 1] isa SecondQuantizedAlgebra.Coeff

        pp = Dict(Δ => 0.0, κ1 => κ_, ϕ => phase)
        ρp = steady(Gphase, pp)
        Rp = frequency_response(Gphase, b, ρp; parameter = pp)
        @test eltype(Rp.scattering) === ComplexF64
        @test scattering_parameter(Rp, 0.0) ≈ -exp(im * phase) atol = 1e-8

        # QIO beam splitters also commonly use ordinary real symbolic sin/cos entries.
        angle = 0.23
        Grotation =
            SLH([cos(ϕ) -sin(ϕ); sin(ϕ) cos(ϕ)], [√(κ1) * a, √(κ2) * a], -Δ * a' * a)
        pr = Dict(Δ => 0.0, κ1 => 0.4, κ2 => 0.6, ϕ => angle)
        ρr = steady(Grotation, pr)
        Rr = frequency_response(Grotation, b, ρr; parameter = pr)
        @test Matrix(Rr.scattering) ≈
              ComplexF64[cos(angle) -sin(angle); sin(angle) cos(angle)]
    end

    @testset "internal loss reduces observable reflection" begin
        κe = 0.3
        κi = 0.7
        p = Dict(Δ => 0.0, κ1 => κe, κ2 => κi)
        ρ = steady(Glin, p)
        R = frequency_response(Glin, b, ρ; parameter = p)
        S11 = scattering_parameter(R, 0.0; in_port = 1, out_port = 1)
        @test S11 ≈ 1 - 2κe / (κe + κi) atol = 1e-8
        @test abs(S11) < 1
    end

    @testset "degenerate parametric amplifier: normal and anomalous scattering" begin
        κ_ = 1.0
        G_ = 0.2
        Gdpa = SLH(1, √(κ1) * a, (Gp / 2) * (a' * a' + a * a))
        p = Dict(κ1 => κ_, Gp => G_)
        ρ = steady(Gdpa, p)
        R = frequency_response(Gdpa, b, ρ; parameter = p)

        analytic(w) = begin
            z = κ_ / 2 - im * w
            D = z^2 - G_^2
            (normal = 1 - κ_ * z / D, anomalous = im * κ_ * G_ / D)
        end

        for w in (-0.7, 0.0, 0.4)
            expected = analytic(w)
            response = scattering_response(R, w)
            @test response.normal[1, 1] ≈ expected.normal rtol = 2e-4 atol = 2e-6
            @test response.anomalous[1, 1] ≈ expected.anomalous rtol = 2e-4 atol = 2e-6
            @test scattering_parameter(R, w; component = :normal) ≈ expected.normal rtol =
                2e-4
            @test scattering_parameter(R, w; component = :anomalous) ≈ expected.anomalous rtol =
                2e-4
            @test abs2(response.normal[1, 1]) - abs2(response.anomalous[1, 1]) ≈ 1 rtol =
                5e-4
        end

        ω = [-1.0, 0.0, 1.0]
        sqz = quadrature_spectrum(R, ω; angle = π / 4)
        anti = quadrature_spectrum(R, ω; angle = 3π / 4)
        analytic_sqz(w) = 1 - 2κ_ * G_ / ((κ_ / 2 + G_)^2 + w^2)
        analytic_anti(w) = 1 + 2κ_ * G_ / ((κ_ / 2 - G_)^2 + w^2)
        @test all(isapprox.(sqz, analytic_sqz.(ω); rtol = 5e-4, atol = 1e-5))
        @test all(isapprox.(anti, analytic_anti.(ω); rtol = 5e-4, atol = 1e-5))
        @test minimum(sqz) < 1
        @test maximum(anti) > 1

        emission = emission_spectrum(R, ω)
        analytic_emission(w) = abs2(analytic(w).anomalous)
        @test all(isapprox.(emission, analytic_emission.(ω); rtol = 5e-4, atol = 1e-5))
        @test all(emission .>= -1e-10)
    end

    @testset "nonlinear Kerr parametric oscillator" begin
        bk = FockBasis(10)
        Hkpo = -0.15 * a' * a + 0.04 * a' * a' * a * a + 0.06 * (a' * a' + a * a)
        Gkpo = SLH(1, a, Hkpo)
        Hn, Jn = to_numeric(Gkpo, bk)
        ρ = steadystate.eigenvector(Hn, collect(Jn))
        R = frequency_response(Gkpo, bk, ρ)
        ω = [-0.4, 0.0, 0.4]

        response = scattering_response(R, ω)
        @test all(isfinite, abs.(response.normal))
        @test all(isfinite, abs.(response.anomalous))
        @test scattering_parameter(R, ω) ≈ collect(response.normal[1, 1, :])
        @test scattering_parameter(Gkpo, bk, ρ, ω) ≈ collect(response.normal[1, 1, :])

        χ = susceptibility(R, a, a', ω)
        @test all(isfinite, abs.(χ))

        emission = emission_spectrum(R, ω)
        @test all(isfinite, emission)
        @test all(emission .>= -1e-9)
    end

    @testset "API contracts and one-shot wrappers" begin
        @test_throws ArgumentError DenseHessenberg(; deflation = 0.0)
        @test_throws ArgumentError DenseHessenberg(; atol = 0.0)

        p = Dict(Δ => 0.0, κ1 => 0.4, κ2 => 0.6)
        ρ = steady(Glin, p)
        R = frequency_response(Glin, b, ρ; parameter = p)
        ω = [-0.3, 0.0, 0.3]

        @test_throws ArgumentError output_field(Glin, 1; input = 1.0)
        @test_throws ArgumentError output_field(R, 1; input = 1.0)
        @test_throws ArgumentError output_field(Glin, 1; input = :invalid)
        @test_throws ArgumentError output_field(R, 1; input = :invalid)
        @test_throws DimensionMismatch output_field(Glin, 1; input = [1.0])
        @test_throws DimensionMismatch output_field(R, 1; input = [1.0])

        @test_throws ArgumentError scattering_response(R, Float64[])
        @test_throws BoundsError scattering_parameter(R, 0.0; in_port = 0)
        @test_throws BoundsError scattering_parameter(R, 0.0; out_port = 3)
        @test_throws ArgumentError scattering_parameter(R, 0.0; component = :invalid)
        @test_throws BoundsError emission_spectrum(R, 0.0; port = 0)
        @test_throws BoundsError quadrature_spectrum(R, 0.0; port = 3)

        @test scattering_response(Glin, b, ρ, ω; parameter = p).normal ≈
              scattering_response(R, ω).normal
        @test susceptibility(Glin, b, ρ, a, a', ω; parameter = p) ≈
              susceptibility(R, a, a', ω)
        @test emission_spectrum(Glin, b, ρ, ω; parameter = p) ≈ emission_spectrum(R, ω)
        @test quadrature_spectrum(Glin, b, ρ, ω; parameter = p, angle = 0.2) ≈
              quadrature_spectrum(R, ω; angle = 0.2)

        @test_throws ArgumentError frequency_response(Glin, b, 2 * ρ; parameter = p)
        excited = dm(fockstate(b, 1))
        @test_throws ArgumentError frequency_response(Glin, b, excited; parameter = p)
        unchecked = frequency_response(
            Glin,
            b,
            excited;
            parameter = p,
            solver = DenseHessenberg(; check_stationary = false),
        )
        @test unchecked isa FrequencyResponse

        mismatched_basis = QuantumOpticsBase.GenericBasis(length(b))
        mismatched_ρ = QuantumOpticsBase.Operator(mismatched_basis, Matrix(ρ.data))
        @test_throws ArgumentError frequency_response(Glin, b, mismatched_ρ; parameter = p)

        a_num = to_numeric(a, b)
        mismatched_a = QuantumOpticsBase.Operator(mismatched_basis, Matrix(a_num.data))
        @test_throws ArgumentError susceptibility(R, mismatched_a, a_num, 0.0)
        @test_throws ArgumentError susceptibility(R, a_num, mismatched_a, 0.0)

        Gone = SLH(1, √(κ1) * a, -Δ * a' * a)
        p1 = Dict(Δ => 0.0, κ1 => 1.0)
        ρ1 = steady(Gone, p1)
        R1 = frequency_response(Gone, b, ρ1; parameter = p1)
        L1 = to_numeric(√(κ1) * a, b; parameter = p1)
        expected = L1 + 2.0 * identityoperator(b)
        @test Matrix(output_field(R1, 1; input = 2.0).data) ≈ Matrix(expected.data)
        @test Matrix(
            to_numeric(output_field(Gone, 1; input = 2.0), b; parameter = p1).data,
        ) ≈ Matrix(expected.data)
    end
end

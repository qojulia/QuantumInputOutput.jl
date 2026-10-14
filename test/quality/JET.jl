using JET
using QuantumInputOutput
using QuantumOptics
using QuantumOpticsBase: dagger
using SecondQuantizedAlgebra
using Test

@testset "JET" begin
    rep = report_package(
        QuantumInputOutput;
        target_modules = (QuantumInputOutput,),
        ignore_missing_comparison = true,
    )
    @show rep
    @test isempty(JET.get_reports(rep))

    @testset "prepared response hot paths" begin
        h = FockSpace(:jet_response)
        a = Destroy(h, :a)
        b = FockBasis(4)
        G = SLH(1, a, 0.2 * a' * a)

        H, J = to_numeric(G, b)
        ρ = steadystate.eigenvector(H, collect(J))
        R = frequency_response(G, b, ρ)
        a_num = to_numeric(a, b)
        adag_num = dagger(a_num)
        omega = [-0.1, 0.0, 0.1]

        JET.@test_opt target_modules = (QuantumInputOutput,) scattering_parameter(R, 0.1)
        JET.@test_opt target_modules = (QuantumInputOutput,) scattering_response(R, 0.1)
        JET.@test_opt target_modules = (QuantumInputOutput,) susceptibility(
            R,
            a_num,
            adag_num,
            0.1,
        )
        JET.@test_opt target_modules = (QuantumInputOutput,) emission_spectrum(R, 0.1)
        JET.@test_opt target_modules = (QuantumInputOutput,) quadrature_spectrum(
            R,
            0.1;
            angle = 0.2,
        )

        JET.@test_opt target_modules = (QuantumInputOutput,) scattering_parameter(R, omega)
        JET.@test_opt target_modules = (QuantumInputOutput,) scattering_response(R, omega)
        JET.@test_opt target_modules = (QuantumInputOutput,) susceptibility(
            R,
            a_num,
            adag_num,
            omega,
        )
        JET.@test_opt target_modules = (QuantumInputOutput,) emission_spectrum(R, omega)
        JET.@test_opt target_modules = (QuantumInputOutput,) quadrature_spectrum(
            R,
            omega;
            angle = 0.2,
        )
    end
end

using JET
using QuantumInputOutput
using Test

@testset "JET" begin
    rep = report_package(
        QuantumInputOutput;
        target_modules = (QuantumInputOutput,),
        ignore_missing_comparison = true,
    )
    @show rep
    @test isempty(JET.get_reports(rep))
end

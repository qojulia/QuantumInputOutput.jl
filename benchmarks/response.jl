"""
Stationary frequency-response benchmarks.
"""
function benchmark_response!(SUITE)
    SUITE["Response"] = BenchmarkGroup()

    h = FockSpace(:response_bench)
    a = Destroy(h, :a)
    @variables κ::Real G::Real

    network = SLH(1, √(κ) * a, (G / 2) * (a' * a' + a * a))
    basis = FockBasis(12)
    parameter = Dict(κ => 1.0, G => 0.2)
    H, J = to_numeric(network, basis; parameter = parameter)
    ρ = steadystate.eigenvector(H, collect(J))
    response = frequency_response(network, basis, ρ; parameter = parameter)
    omega = collect(range(-2.0, 2.0; length = 201))

    SUITE["Response"]["prepare"] =
        @benchmarkable frequency_response($network, $basis, $ρ; parameter = $parameter)
    SUITE["Response"]["S11 sweep prepared"] =
        @benchmarkable scattering_parameter($response, $omega)
    SUITE["Response"]["full Bogoliubov sweep prepared"] =
        @benchmarkable scattering_response($response, $omega)
    SUITE["Response"]["emission sweep prepared"] =
        @benchmarkable emission_spectrum($response, $omega)
    SUITE["Response"]["quadrature sweep prepared"] =
        @benchmarkable quadrature_spectrum($response, $omega; angle = π / 4)

    return nothing
end

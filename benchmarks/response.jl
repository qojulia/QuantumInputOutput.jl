struct _LegacyResponseResolvent{BT,FT}
    basis::BT
    factorization::FT
    d::Int
end

function _legacy_response_resolvent(H, J, ρ)
    basis = H.basis_l
    d = length(basis)
    ℒ = QuantumOpticsBase.liouvillian(H, collect(J))
    vecI = reshape(Matrix(QuantumOpticsBase.identityoperator(basis).data), d * d)
    vecρ = reshape(Matrix(ρ.data), d * d)
    factorization = LinearAlgebra.hessenberg(Matrix(ℒ.data) + vecρ * vecI')
    return _LegacyResponseResolvent(basis, factorization, d)
end

function _legacy_response_solve(R::_LegacyResponseResolvent, μ::Number, rhs)
    x = (R.factorization + μ * LinearAlgebra.I) \ reshape(Matrix(rhs.data), R.d * R.d)
    return QuantumOpticsBase.Operator(R.basis, R.basis, reshape(x, R.d, R.d))
end

function _legacy_scattering_sweep(R::_LegacyResponseResolvent, Lout, rhs, direct, omega)
    return [
        direct + LinearAlgebra.tr(Lout.data * _legacy_response_solve(R, im * ω, rhs).data)
        for ω in omega
    ]
end

"""
Stationary frequency-response benchmarks.

The `legacy operator sweep` reproduces the hot loop used by the original #21
prototype after its Hessenberg factorization has already been prepared. It
therefore isolates the cost removed by the new Liouville-vector cache: temporary
`Operator` reconstruction and matrix-product trace contractions at every
frequency.
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

    legacy = _legacy_response_resolvent(H, J, ρ)
    Lout = J[1]
    rhs = dagger(J[1]) * ρ - ρ * dagger(J[1])
    legacy_reference = _legacy_scattering_sweep(legacy, Lout, rhs, 1.0 + 0im, omega)
    prepared_reference = scattering_parameter(response, omega)
    @assert maximum(abs.(legacy_reference .- prepared_reference)) < 1e-10

    SUITE["Response"]["prepare"] =
        @benchmarkable frequency_response($network, $basis, $ρ; parameter = $parameter)
    SUITE["Response"]["S11 single prepared"] =
        @benchmarkable scattering_parameter($response, 0.1)
    SUITE["Response"]["S11 sweep prepared"] =
        @benchmarkable scattering_parameter($response, $omega)
    SUITE["Response"]["S11 sweep legacy operator"] =
        @benchmarkable _legacy_scattering_sweep($legacy, $Lout, $rhs, $(1.0 + 0im), $omega)
    SUITE["Response"]["S11 sweep one-shot"] = @benchmarkable scattering_parameter(
        $network,
        $basis,
        $ρ,
        $omega;
        parameter = $parameter,
    )
    SUITE["Response"]["full Bogoliubov sweep prepared"] =
        @benchmarkable scattering_response($response, $omega)
    SUITE["Response"]["emission sweep prepared"] =
        @benchmarkable emission_spectrum($response, $omega)
    SUITE["Response"]["quadrature sweep prepared"] =
        @benchmarkable quadrature_spectrum($response, $omega; angle = π / 4)

    return nothing
end

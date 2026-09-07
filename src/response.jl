################################################################
### Frequency-domain response of stationary SLH networks      ###
################################################################

abstract type AbstractResponseSolver end

"""
    DenseHessenberg(; deflation=1.0, check_stationary=true, atol=1e-10)

Dense shifted-resolvent backend for stationary response calculations.

The Liouvillian is converted to a dense matrix once, its unique stationary zero
mode is deflated with `ρ_ss`, and a Hessenberg factorisation is reused for all
frequency shifts. This is intended for small and medium Hilbert spaces where a
dense Liouville-space representation is appropriate.

`check_stationary=true` verifies that `ρ_ss` has unit trace and is stationary to
`atol`. The backend assumes that the stationary state is unique; genuine
additional zero modes or undamped dynamical modes are not regularised by the
rank-one deflation.
"""
struct DenseHessenberg <: AbstractResponseSolver
    deflation::Float64
    check_stationary::Bool
    atol::Float64
end

function DenseHessenberg(;
    deflation::Real = 1.0,
    check_stationary::Bool = true,
    atol::Real = 1e-10,
)
    iszero(deflation) && throw(ArgumentError("`deflation` must be nonzero"))
    atol > 0 || throw(ArgumentError("`atol` must be positive"))
    return DenseHessenberg(Float64(deflation), check_stationary, Float64(atol))
end

struct HessenbergResolvent{FT}
    factorization::FT
end

struct PortResponseCache{SRT,OWT,ERT,FWT,FDWT}
    scattering_rhs::SRT
    output_weights::OWT
    emission_rhs::ERT
    fluctuation_weights::FWT
    fluctuation_dagger_weights::FDWT
end

"""
    FrequencyResponse

Prepared stationary response problem for an SLH network.

Create one with [`frequency_response`](@ref), then reuse it for
[`scattering_response`](@ref), [`scattering_parameter`](@ref),
[`susceptibility`](@ref), [`emission_spectrum`](@ref), and
[`quadrature_spectrum`](@ref). The translated Hamiltonian, jump operators,
steady state, port-response vectors, and shifted-resolvent backend are prepared
once and reused by all subsequent observables.
"""
struct FrequencyResponse{GT,BT,PT,OPT,ST,HT,JT,RT,RDT,BET,CT}
    network::GT
    basis::BT
    parameter::PT
    operators::OPT
    scattering::ST
    hamiltonian::HT
    jump_operators::JT
    steady_state::RT
    rho_data::RDT
    backend::BET
    cache::CT
end

"""
    ScatteringResponse

Normal and anomalous weak-probe scattering response.

For a scalar frequency, `normal` and `anomalous` are port-by-port matrices. For
a vector of frequencies, the third array dimension indexes frequency. The
normal block maps input annihilation amplitudes to output annihilation
amplitudes; the anomalous block maps conjugate input amplitudes to output
annihilation amplitudes. The anomalous block vanishes for passive
number-conserving systems but is generally nonzero for parametrically driven
systems.
"""
struct ScatteringResponse{WT,NT,AT}
    omega::WT
    normal::NT
    anomalous::AT
end

# ──────────────────────────────────────────────
# Parameter/scattering preparation
# ──────────────────────────────────────────────

_response_parameter_value(
    x::Number,
    parameter,
    ρ_ss::QuantumOpticsBase.AbstractOperator,
)::ComplexF64 = ComplexF64(x)

function _response_parameter_value(
    x::Complex{Symbolics.Num},
    parameter,
    ρ_ss::QuantumOpticsBase.AbstractOperator,
)::ComplexF64
    return SQA.numeric_average(substitute(x, parameter), ρ_ss)
end

function _response_parameter_value(
    x::Symbolics.Num,
    parameter,
    ρ_ss::QuantumOpticsBase.AbstractOperator,
)::ComplexF64
    x_ = substitute(x, parameter)
    return SQA.numeric_average(complex(x_, zero(x_)), ρ_ss)
end

function _response_parameter_value(
    x::SQA.Coeff,
    parameter,
    ρ_ss::QuantumOpticsBase.AbstractOperator,
)::ComplexF64
    x_ = SQA.substitute(x, parameter)
    return SQA.numeric_average(SQA.to_num(x_), ρ_ss)
end

function _response_parameter_value(
    x::BasicSymbolic,
    parameter,
    ρ_ss::QuantumOpticsBase.AbstractOperator,
)::ComplexF64
    return SQA.numeric_average(substitute(x, parameter), ρ_ss)
end

function _response_parameter_value(x, parameter, ρ_ss)
    throw(
        ArgumentError(
            "frequency-domain response requires a time-independent numeric scattering matrix; " *
            "could not evaluate entry of type $(typeof(x))",
        ),
    )
end

function _numeric_scattering(G::SLH{N}, parameter, ρ_ss) where {N}
    S = scattering(G)
    values = ntuple(k -> _response_parameter_value(S[k], parameter, ρ_ss), Val(N * N))
    return SMatrix{N,N,ComplexF64}(values)
end

# ──────────────────────────────────────────────
# Dense Liouville-space helpers
# ──────────────────────────────────────────────

_operator_matrix(A::QuantumOpticsBase.AbstractOperator) = Matrix(A.data)

function _trace_weights(A::AbstractMatrix)
    return vec(Matrix(transpose(A)))
end

function _contract(weights::AbstractVector, x::AbstractVector)
    T = promote_type(eltype(weights), eltype(x))
    value = zero(T)
    @inbounds for i in eachindex(weights, x)
        value += weights[i] * x[i]
    end
    return value
end

function _commutator_vector(A::AbstractMatrix, ρ::AbstractMatrix)
    return vec(A * ρ - ρ * A)
end

_left_product_vector(A::AbstractMatrix, ρ::AbstractMatrix) = vec(A * ρ)

function _fluctuation_matrix(A::AbstractMatrix, ρ::AbstractMatrix)
    mean_A = _contract(_trace_weights(A), vec(ρ))
    δA = copy(A)
    @inbounds for i in axes(δA, 1)
        δA[i, i] -= mean_A
    end
    return δA
end

function _trace_vector(::Type{T}, d::Int) where {T}
    vecI = zeros(T, d * d)
    @inbounds for i = 1:d
        vecI[i+(i-1)*d] = one(T)
    end
    return vecI
end

function _prepare_backend(
    H::QuantumOpticsBase.AbstractOperator,
    J,
    ρ_ss::QuantumOpticsBase.AbstractOperator,
    solver::DenseHessenberg,
)
    ℒ = QuantumOpticsBase.liouvillian(H, collect(J))
    Lmat = Matrix(ℒ.data)
    ρmat = _operator_matrix(ρ_ss)
    d = size(ρmat, 1)
    size(ρmat, 2) == d || throw(ArgumentError("`ρ_ss` must be square"))
    size(Lmat) == (d * d, d * d) ||
        throw(DimensionMismatch("Liouvillian and steady-state dimensions do not agree"))

    vecρ = vec(copy(ρmat))
    vecI = _trace_vector(eltype(Lmat), d)

    if solver.check_stationary
        trace_error = abs(_contract(vecI, vecρ) - one(eltype(vecρ)))
        trace_error <= solver.atol ||
            throw(ArgumentError("`ρ_ss` is not normalized: |tr(ρ_ss)-1| = $trace_error"))

        residual = LinearAlgebra.norm(Lmat * vecρ)
        residual <= solver.atol ||
            throw(ArgumentError("`ρ_ss` is not stationary: ||ℒρ_ss|| = $residual"))
    end

    # Rank-one deflation ℒ -> ℒ + η |ρ_ss><1|. Mutate the dense matrix in place
    # rather than materialising an additional d² × d² outer-product matrix.
    η = solver.deflation
    @inbounds for j in eachindex(vecI)
        c = vecI[j]
        iszero(c) && continue
        for i in eachindex(vecρ)
            Lmat[i, j] += η * vecρ[i] * c
        end
    end

    F = LinearAlgebra.hessenberg!(Lmat)
    return HessenbergResolvent(F)
end

function _shifted_solve(R::HessenbergResolvent, μ::Number, rhs)
    return (R.factorization + μ * LinearAlgebra.I) \ rhs
end

# ──────────────────────────────────────────────
# Port caches
# ──────────────────────────────────────────────

function _input_coupling_matrix(S, Jm, input_port::Int)
    B = conj(S[1, input_port]) * Jm[1]
    @inbounds for k = 2:length(Jm)
        B .+= conj(S[k, input_port]) .* Jm[k]
    end
    return B
end

function _prepare_port_cache(S, J, ρmat)
    nports = length(J)
    nports > 0 || throw(ArgumentError("frequency response requires at least one port"))
    Jm = [_operator_matrix(Jk) for Jk in J]

    output_weights = hcat((_trace_weights(A) for A in Jm)...)
    input_couplings = [_input_coupling_matrix(S, Jm, j) for j = 1:nports]
    normal_rhs = hcat((_commutator_vector(adjoint(B), ρmat) for B in input_couplings)...)
    anomalous_rhs = hcat((-_commutator_vector(B, ρmat) for B in input_couplings)...)
    scattering_rhs = hcat(normal_rhs, anomalous_rhs)

    δJ = [_fluctuation_matrix(A, ρmat) for A in Jm]
    emission_rhs = hcat((_left_product_vector(A, ρmat) for A in δJ)...)
    fluctuation_weights = hcat((_trace_weights(A) for A in δJ)...)
    fluctuation_dagger_weights = hcat((_trace_weights(adjoint(A)) for A in δJ)...)

    return PortResponseCache(
        scattering_rhs,
        output_weights,
        emission_rhs,
        fluctuation_weights,
        fluctuation_dagger_weights,
    )
end

# ──────────────────────────────────────────────
# Prepared response construction
# ──────────────────────────────────────────────

"""
    frequency_response(G::SLH, basis, ρ_ss; parameter=Dict(), operators=Dict(),
                       solver=DenseHessenberg())

Prepare the stationary frequency-domain response of an SLH network around the
steady state `ρ_ss`.

The network is translated once with [`to_numeric`](@ref), including the numeric
working-point value of its scattering matrix. The response backend and all
port-dependent Liouville-space vectors are then cached for reuse.

This interface currently supports time-independent Hamiltonians, jump operators,
and scattering matrices. The dense default solver assumes a unique stationary
state. `ρ_ss` may be obtained with any steady-state method compatible with the
same numeric model.
"""
function frequency_response(
    G::SLH,
    basis,
    ρ_ss::QuantumOpticsBase.AbstractOperator;
    parameter = Dict(),
    operators = Dict(),
    solver::AbstractResponseSolver = DenseHessenberg(),
)
    parameter_ = copy(parameter)
    operators_ = copy(operators)
    H, J = to_numeric(G, basis; parameter = parameter_, operators = operators_)

    H isa QuantumOpticsBase.AbstractOperator ||
        throw(ArgumentError("frequency response requires a time-independent Hamiltonian"))
    all(Jk -> Jk isa QuantumOpticsBase.AbstractOperator, J) ||
        throw(ArgumentError("frequency response requires time-independent jump operators"))

    S = _numeric_scattering(G, parameter_, ρ_ss)
    ρmat = _operator_matrix(ρ_ss)
    backend = _prepare_backend(H, J, ρ_ss, solver)
    cache = _prepare_port_cache(S, J, ρmat)

    return FrequencyResponse(
        G,
        basis,
        parameter_,
        operators_,
        S,
        H,
        J,
        ρ_ss,
        ρmat,
        backend,
        cache,
    )
end

function _prepare_backend(H, J, ρ_ss, solver::AbstractResponseSolver)
    throw(ArgumentError("unsupported frequency-response solver $(typeof(solver))"))
end

# ──────────────────────────────────────────────
# Output field
# ──────────────────────────────────────────────

"""
    output_field(G::SLH, port; input=nothing)

Return the system-operator contribution to the travelling output field at
`port`. QuantumInputOutput uses the standard SLH convention

```math
b_{\\mathrm{out}} = S b_{\\mathrm{in}} + L.
```

With `input=nothing`, this returns `jump_operator(G)[port]`. Passing a scalar
coherent input for a one-port network or a vector of coherent input amplitudes
adds the classical directly scattered contribution `(S*input)[port]`.
"""
function output_field(G::SLH, port::Integer; input = nothing)
    L = jump_operator(G)[port]
    input === nothing && return L
    S = scattering(G)
    nports = size(S, 2)
    if input isa Number
        nports == 1 ||
            throw(ArgumentError("scalar `input` is only valid for a one-port network"))
        return S[port, 1] * input + L
    end
    input isa AbstractVector || throw(ArgumentError("`input` must be a scalar or vector"))
    length(input) == nports ||
        throw(DimensionMismatch("input amplitude vector has the wrong length"))
    offset = sum(S[port, j] * input[j] for j in axes(S, 2))
    return offset + L
end

function output_field(R::FrequencyResponse, port::Integer; input = nothing)
    _check_port(R, port)
    L = R.jump_operators[port]
    input === nothing && return L
    nports = length(R.jump_operators)
    if input isa Number
        nports == 1 ||
            throw(ArgumentError("scalar `input` is only valid for a one-port network"))
        offset = R.scattering[port, 1] * input
    else
        input isa AbstractVector ||
            throw(ArgumentError("`input` must be a scalar or vector"))
        length(input) == nports ||
            throw(DimensionMismatch("input amplitude vector has the wrong length"))
        offset = sum(R.scattering[port, j] * input[j] for j in axes(R.scattering, 2))
    end
    return L + offset * QuantumOpticsBase.identityoperator(R.basis)
end

# ──────────────────────────────────────────────
# Susceptibility
# ──────────────────────────────────────────────

function _numeric_response_operator(
    R::FrequencyResponse,
    A::QuantumOpticsBase.AbstractOperator,
)
    return A
end

function _numeric_response_operator(R::FrequencyResponse, A)
    A_ = to_numeric(A, R.basis; parameter = R.parameter, operators = R.operators)
    A_ isa QuantumOpticsBase.AbstractOperator ||
        throw(ArgumentError("response observables must be time-independent operators"))
    return A_
end

function _susceptibility_data(R::FrequencyResponse, A, B)
    A_ = _numeric_response_operator(R, A)
    B_ = _numeric_response_operator(R, B)
    weights = _trace_weights(_operator_matrix(A_))
    rhs = _commutator_vector(_operator_matrix(B_), R.rho_data)
    return weights, rhs
end

function _susceptibility_at(R::FrequencyResponse, weights, rhs, ω::Real)
    x = _shifted_solve(R.backend, im * ω, rhs)
    return im * _contract(weights, x)
end

"""
    susceptibility(R::FrequencyResponse, A, B, omega)
    susceptibility(G::SLH, basis, ρ_ss, A, B, omega; kwargs...)

Kubo susceptibility of observable `A` to a weak perturbation coupled through
`B`,

```math
\\chi_{AB}(\\omega) = i\\,\\mathrm{Tr}\\left[
A(\\mathcal L+i\\omega)^{-1}[B,\\rho_{ss}]\\right].
```

`A` and `B` may be numeric operators or symbolic operators compatible with the
parameter/operator context stored in the prepared [`FrequencyResponse`](@ref).
`omega` may be a scalar or an `AbstractVector`.
"""
function susceptibility(R::FrequencyResponse, A, B, omega::Real)
    weights, rhs = _susceptibility_data(R, A, B)
    return _susceptibility_at(R, weights, rhs, omega)
end

function susceptibility(R::FrequencyResponse, A, B, omega::AbstractVector)
    weights, rhs = _susceptibility_data(R, A, B)
    return [_susceptibility_at(R, weights, rhs, ω) for ω in omega]
end

function susceptibility(G::SLH, basis, ρ_ss, A, B, omega; kwargs...)
    R = frequency_response(G, basis, ρ_ss; kwargs...)
    return susceptibility(R, A, B, omega)
end

# ──────────────────────────────────────────────
# Normal/anomalous multiport scattering
# ──────────────────────────────────────────────

_nports(R::FrequencyResponse) = length(R.jump_operators)

function _check_port(R::FrequencyResponse, port::Integer)
    1 <= port <= _nports(R) || throw(BoundsError(R.jump_operators, port))
    return nothing
end

function _scattering_at(R::FrequencyResponse, ω::Real)
    nports = _nports(R)
    X = _shifted_solve(R.backend, im * ω, R.cache.scattering_rhs)
    dynamic = transpose(R.cache.output_weights) * X
    normal = Matrix(R.scattering) + dynamic[:, 1:nports]
    anomalous = dynamic[:, (nports+1):(2*nports)]
    return normal, anomalous
end

"""
    scattering_response(R::FrequencyResponse, omega)
    scattering_response(G::SLH, basis, ρ_ss, omega; kwargs...)

Return the complete weak-probe multiport scattering response.

For parametrically driven systems the response is Bogoliubov/Nambu-like: the
`normal` block maps input annihilation amplitudes to output annihilation
amplitudes and `anomalous` maps conjugate input amplitudes to output
annihilation amplitudes. For a scalar frequency both fields are matrices. For a
frequency vector they are `nport × nport × nfrequency` arrays.

The drive seen by the internal system is `(S' * L)_j`, so nontrivial static SLH
scattering matrices are included in both the direct and dynamical response.
"""
function scattering_response(R::FrequencyResponse, omega::Real)
    normal, anomalous = _scattering_at(R, omega)
    return ScatteringResponse(omega, normal, anomalous)
end

function scattering_response(R::FrequencyResponse, omega::AbstractVector)
    isempty(omega) && throw(ArgumentError("`omega` must not be empty"))
    first_normal, first_anomalous = _scattering_at(R, first(omega))
    nports = _nports(R)
    nfreq = length(omega)
    T = promote_type(eltype(first_normal), eltype(first_anomalous))
    normal = Array{T}(undef, nports, nports, nfreq)
    anomalous = Array{T}(undef, nports, nports, nfreq)
    normal[:, :, 1] = first_normal
    anomalous[:, :, 1] = first_anomalous
    @inbounds for k = 2:nfreq
        normal_k, anomalous_k = _scattering_at(R, omega[k])
        normal[:, :, k] = normal_k
        anomalous[:, :, k] = anomalous_k
    end
    return ScatteringResponse(omega, normal, anomalous)
end

function scattering_response(G::SLH, basis, ρ_ss, omega; kwargs...)
    R = frequency_response(G, basis, ρ_ss; kwargs...)
    return scattering_response(R, omega)
end

function _scattering_parameter_at(
    R::FrequencyResponse,
    ω::Real,
    in_port::Integer,
    out_port::Integer,
    component::Symbol,
)
    _check_port(R, in_port)
    _check_port(R, out_port)
    nports = _nports(R)
    column = if component === :normal
        in_port
    elseif component === :anomalous
        nports + in_port
    else
        throw(ArgumentError("`component` must be `:normal` or `:anomalous`"))
    end
    rhs = view(R.cache.scattering_rhs, :, column)
    x = _shifted_solve(R.backend, im * ω, rhs)
    dynamic = _contract(view(R.cache.output_weights, :, out_port), x)
    return component === :normal ? R.scattering[out_port, in_port] + dynamic : dynamic
end

"""
    scattering_parameter(R::FrequencyResponse, omega; in_port=1, out_port=in_port,
                         component=:normal)

Return one element of the prepared weak-probe scattering response. `component`
may be `:normal` or `:anomalous`. If `out_port` is omitted it defaults to
`in_port`, so the natural one-port operation is reflection.
"""
function scattering_parameter(
    R::FrequencyResponse,
    omega::Real;
    in_port::Integer = 1,
    out_port::Union{Nothing,Integer} = nothing,
    component::Symbol = :normal,
)
    out = isnothing(out_port) ? in_port : out_port
    return _scattering_parameter_at(R, omega, in_port, out, component)
end

function scattering_parameter(
    R::FrequencyResponse,
    omega::AbstractVector;
    in_port::Integer = 1,
    out_port::Union{Nothing,Integer} = nothing,
    component::Symbol = :normal,
)
    out = isnothing(out_port) ? in_port : out_port
    return [_scattering_parameter_at(R, ω, in_port, out, component) for ω in omega]
end

function scattering_parameter(
    G::SLH,
    basis,
    ρ_ss,
    omega;
    parameter = Dict(),
    operators = Dict(),
    solver::AbstractResponseSolver = DenseHessenberg(),
    in_port::Integer = 1,
    out_port::Union{Nothing,Integer} = nothing,
    component::Symbol = :normal,
)
    R = frequency_response(
        G,
        basis,
        ρ_ss;
        parameter = parameter,
        operators = operators,
        solver = solver,
    )
    return scattering_parameter(
        R,
        omega;
        in_port = in_port,
        out_port = out_port,
        component = component,
    )
end

# ──────────────────────────────────────────────
# Output spectra
# ──────────────────────────────────────────────

function _emission_at(R::FrequencyResponse, ω::Real, port::Integer)
    _check_port(R, port)
    rhs = view(R.cache.emission_rhs, :, port)
    weights = view(R.cache.fluctuation_dagger_weights, :, port)
    x = _shifted_solve(R.backend, -im * ω, rhs)
    return -2 * real(_contract(weights, x))
end

"""
    emission_spectrum(R::FrequencyResponse, omega; port=1)

Normally ordered connected output-emission spectrum at `port`. The vacuum floor
is zero. The coherent stationary mean is removed before the resolvent solve;
its elastic contribution is a delta peak and is intentionally not folded into
the smooth spectrum.
"""
emission_spectrum(R::FrequencyResponse, omega::Real; port::Integer = 1) =
    _emission_at(R, omega, port)

function emission_spectrum(R::FrequencyResponse, omega::AbstractVector; port::Integer = 1)
    return [_emission_at(R, ω, port) for ω in omega]
end

function emission_spectrum(G::SLH, basis, ρ_ss, omega; port::Integer = 1, kwargs...)
    R = frequency_response(G, basis, ρ_ss; kwargs...)
    return emission_spectrum(R, omega; port = port)
end

function _quadrature_at(R::FrequencyResponse, ω::Real, port::Integer, angle::Real)
    _check_port(R, port)
    rhs = view(R.cache.emission_rhs, :, port)
    wL = view(R.cache.fluctuation_weights, :, port)
    wLd = view(R.cache.fluctuation_dagger_weights, :, port)

    xplus = _shifted_solve(R.backend, im * ω, rhs)
    xminus = iszero(ω) ? xplus : _shifted_solve(R.backend, -im * ω, rhs)

    corr_ll_plus = -_contract(wL, xplus)
    corr_ll_minus = -_contract(wL, xminus)
    corr_ldl_plus = -_contract(wLd, xplus)
    corr_ldl_minus = -_contract(wLd, xminus)

    return 1 +
           2 * real(
        exp(-2im * angle) * (corr_ll_plus + corr_ll_minus) + corr_ldl_plus + corr_ldl_minus,
    )
end

"""
    quadrature_spectrum(R::FrequencyResponse, omega; port=1, angle=0)

Vacuum-normalised homodyne spectrum of the connected output quadrature

```math
X_\\theta = (e^{-i\\theta}b_{out} + e^{i\\theta}b_{out}^\\dagger)/\\sqrt{2}.
```

The shot-noise floor is one. A squeezed quadrature therefore falls below one.
Only two shifted resolvent solves are required per nonzero frequency because the
normal and anomalous contractions share the same right-hand side.
"""
function quadrature_spectrum(
    R::FrequencyResponse,
    omega::Real;
    port::Integer = 1,
    angle::Real = 0,
)
    return _quadrature_at(R, omega, port, angle)
end

function quadrature_spectrum(
    R::FrequencyResponse,
    omega::AbstractVector;
    port::Integer = 1,
    angle::Real = 0,
)
    return [_quadrature_at(R, ω, port, angle) for ω in omega]
end

function quadrature_spectrum(
    G::SLH,
    basis,
    ρ_ss,
    omega;
    port::Integer = 1,
    angle::Real = 0,
    kwargs...,
)
    R = frequency_response(G, basis, ρ_ss; kwargs...)
    return quadrature_spectrum(R, omega; port = port, angle = angle)
end

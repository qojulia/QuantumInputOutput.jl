######################################
### functions for virtual cavities ###
######################################

const _tol_div = 1e-10
const _extrapolate = ExtrapolationType.Extension
const _ϵ = 1e-10

"""
    PulseCoupling(f)

Callable wrapper used for time-dependent pulse couplings. A `PulseCoupling` has exactly one
call method, `g(t::Real)`, and can also be broadcast as `g.(T)`.

The wrapper is useful as a `time_parameter` value for [`to_numeric`](@ref):
SecondQuantizedAlgebra requires such values to have an unambiguous single-argument call
interface, while interpolation objects may expose additional call arities.

`f` may be an analytic closure, an interpolation, or any evaluator accepting one real time.
"""
struct PulseCoupling{F}
    f::F
end
(g::PulseCoupling)(t::Real) = g.f(t)

_interpolate_mode(mode::AbstractVector, T::AbstractVector) =
    LinearInterpolation(mode, T; extrapolation = _extrapolate)
_interpolate_mode(mode, T::AbstractVector) = mode

# ──────────────────────────────────────────────
# Gaussian pulse type
# ──────────────────────────────────────────────

"""
    Gaussian(τ, σ; δ=0)

Describe the normalized Gaussian temporal mode

```math
u(t)=\frac{1}{\sqrt{\sigma}\,\pi^{1/4}}
\exp\!\left[-\frac{(t-\tau)^2}{2\sigma^2}\right]e^{i\delta t}.
```

`τ` is the pulse center, `σ` its width, and `δ` an optional detuning. The descriptor is
used by [`coupling_input`](@ref) and [`coupling_output`](@ref) to evaluate the mode and its
accumulated norm analytically rather than from sampled quadrature. The analytic cumulative
integral uses `t = 0` as the initial time.
"""
struct Gaussian{T}
    τ::T
    σ::T
    δ::T
end
Gaussian(τ, σ; δ = zero(τ)) = Gaussian(promote(τ, σ, δ)...)

function _gaussian_mode(g::Gaussian)
    τ, σ, δ = g.τ, g.σ, g.δ
    mode = if δ == 0
        t -> 1 / (√(σ) * π^(1 / 4)) * exp(-0.5 * (t - τ)^2 / σ^2)
    else
        t -> 1 / (√(σ) * π^(1 / 4)) * exp(-0.5 * (t - τ)^2 / σ^2) * exp(1im * δ * t)
    end
    ∫mode2 = t -> 0.5 * (erf((t - τ) / σ) + erf(τ / σ))
    return mode, ∫mode2
end

# ──────────────────────────────────────────────
# Shared coupling core
# ──────────────────────────────────────────────

function _coupling_from_mode(mode::Vector, T::Vector, denominator)
    nt = length(T)
    mode_sq = Vector{Float64}(undef, nt)
    map!(abs2, mode_sq, mode)
    ∫mode2 = cumul_integrate(T, mode_sq)
    coupling = zeros(ComplexF64, nt)
    @inbounds for i = 1:nt
        denom = denominator(∫mode2[i])
        if sqrt(abs(denom)) > _tol_div
            coupling[i] = mode[i]' / sqrt(denom)
        end
    end
    return PulseCoupling(LinearInterpolation(coupling, T; extrapolation = _extrapolate))
end

# ──────────────────────────────────────────────
# coupling_input / coupling_output
# ──────────────────────────────────────────────

"""
    coupling_input(u, T)

Construct the source-cavity coupling for a normalized input temporal mode `u` on the time
grid `T`. `u` may be a vector sampled on `T`, a function `u(t)`, or a
`LinearInterpolation`.

With

```math
F(t)=\int_{T[1]}^t |u(t')|^2\,dt',
```

the virtual-cavity construction uses

```math
g_u(t)=\frac{u^*(t)}{\sqrt{1-F(t)}}.
```

The sampled implementation evaluates the cumulative integral numerically, regularizes the
endpoint denominator, and returns a [`PulseCoupling`](@ref) backed by linear interpolation.
The coupling can be supplied directly as a `time_parameter` value to [`to_numeric`](@ref).

For a Gaussian mode with an analytic cumulative norm, use `coupling_input(::Gaussian)`.
See [Quantum pulses and input-output theory](@ref) for the virtual-cavity derivation.
"""
coupling_input(u::Vector, T::Vector) = _coupling_from_mode(u, T, x -> abs(1 - x) + _ϵ)
coupling_input(u::Function, T::Vector) = coupling_input(u.(T), T)
coupling_input(u::LinearInterpolation, T::Vector) = coupling_input(u.(T), T)

"""
    coupling_input(g::Gaussian)

Construct the analytic source-cavity coupling for the Gaussian mode `g`. The Gaussian
envelope and its accumulated norm are evaluated in closed form, and the result is returned
as a [`PulseCoupling`](@ref).
"""
function coupling_input(g::Gaussian)
    mode, ∫m2 = _gaussian_mode(g)
    return PulseCoupling(t -> mode(t)' / √(abs(1 - ∫m2(t)) + _ϵ))
end

"""
    coupling_output(v, T)

Construct the receiver-cavity coupling for a normalized output temporal mode `v` on the time
grid `T`. `v` may be a vector sampled on `T`, a function `v(t)`, or a
`LinearInterpolation`.

With

```math
F(t)=\int_{T[1]}^t |v(t')|^2\,dt',
```

the virtual-cavity construction uses

```math
g_v(t)=-\frac{v^*(t)}{\sqrt{F(t)}}.
```

The sampled implementation evaluates the cumulative integral numerically, regularizes the
endpoint denominator, and returns a [`PulseCoupling`](@ref) backed by linear interpolation.
The coupling can be supplied directly as a `time_parameter` value to [`to_numeric`](@ref).

For a Gaussian mode with an analytic cumulative norm, use `coupling_output(::Gaussian)`.
See [Quantum pulses and input-output theory](@ref) for the virtual-cavity derivation.
"""
coupling_output(v::Vector, T::Vector) = _coupling_from_mode(-v, T, x -> x + _ϵ)
coupling_output(v::Function, T::Vector) = coupling_output(v.(T), T)
coupling_output(v::LinearInterpolation, T::Vector) = coupling_output(v.(T), T)

"""
    coupling_output(g::Gaussian)

Construct the analytic receiver-cavity coupling for the Gaussian mode `g`. The Gaussian
envelope and its accumulated norm are evaluated in closed form, and the result is returned
as a [`PulseCoupling`](@ref).
"""
function coupling_output(g::Gaussian)
    mode, ∫m2 = _gaussian_mode(g)
    return PulseCoupling(t -> -mode(t)' / √(∫m2(t) + _ϵ))
end

# ──────────────────────────────────────────────
# effective_output_mode (was v_eff)
# ──────────────────────────────────────────────

"""
    effective_output_mode(v_fcts, gv_fcts, T, i; alg=Tsit5(), kwargs...)
    effective_output_mode(v_fcts, T, i; alg=Tsit5(), kwargs...)

Return the effective temporal mode seen by output cavity `i` when several receiver cavities
are cascaded after the physical system. The preceding receiver cavities distort the traveling
field, so the nominal mode `v_fcts[i]` must be corrected before its coupling is constructed.

The output modes must be ordered starting with the first output cavity after the system, and
`i > 1`. `v_fcts` may contain callables or sampled mode vectors. When `gv_fcts` is omitted,
[`coupling_output`](@ref) is used to construct the couplings of the nominal modes. Sampled
coupling data may also be supplied explicitly.

Returns a callable `vᵢ_eff(t)`. `alg` and all remaining keyword arguments are forwarded to
the auxiliary OrdinaryDiffEq solve.

See also [`effective_input_mode`](@ref) and [Pulse modes](@ref).
"""
function effective_output_mode(v_fcts, gv_fcts, T, i; alg = Tsit5(), kwargs...)
    @assert i > 1
    n = i - 1
    # Capture as a tuple for concrete closure types
    gv = ntuple(k -> gv_fcts[k], n)
    v_i = v_fcts[i]
    function multiple_outputs_α!(dα, α, p, t)
        gv_buf = p
        @inbounds for k = 1:n
            gv_buf[k] = gv[k](t)
        end
        vi_t = v_i(t)
        @inbounds for j = 1:n
            coupling_sum = zero(ComplexF64)
            for k = 1:(j-1)
                coupling_sum += gv_buf[k]' * α[k]
            end
            dα[j] = -gv_buf[j] * (vi_t + coupling_sum) - 0.5 * abs2(gv_buf[j]) * α[j]
        end
    end
    u0 = zeros(ComplexF64, n)
    p = Vector{ComplexF64}(undef, n)
    tspan = (T[1], T[end])
    prob = ODEProblem(multiple_outputs_α!, u0, tspan, p)
    sol_α = solve(prob, alg; kwargs...)
    function v_i_eff(t)
        α_t = sol_α(t)
        result = v_i(t)
        @inbounds for k = 1:n
            result += gv[k](t)' * α_t[k]
        end
        return result
    end
    return v_i_eff
end

effective_output_mode(v_fcts, T, i; kwargs...) = effective_output_mode(
    v_fcts,
    [coupling_output(v_, T) for v_ in v_fcts],
    T,
    i;
    kwargs...,
)

function effective_output_mode(
    v_data::AbstractVector{<:AbstractVector},
    gv_data::AbstractVector{<:AbstractVector},
    T::AbstractVector,
    i;
    alg = Tsit5(),
    kwargs...,
)
    v_fcts = [_interpolate_mode(v_, T) for v_ in v_data]
    gv_fcts = [_interpolate_mode(gv_, T) for gv_ in gv_data]
    return effective_output_mode(v_fcts, gv_fcts, T, i; alg, kwargs...)
end

effective_output_mode(
    v_data::AbstractVector{<:AbstractVector},
    T::AbstractVector,
    i;
    alg = Tsit5(),
    kwargs...,
) = effective_output_mode(
    v_data,
    [coupling_output(v_, T).(T) for v_ in v_data],
    T,
    i;
    alg,
    kwargs...,
)

# ──────────────────────────────────────────────
# effective_input_mode (was u_eff)
# ──────────────────────────────────────────────

"""
    effective_input_mode(u_fcts, gu_fcts, T, i; alg=Tsit5(), kwargs...)
    effective_input_mode(u_fcts, T, i; alg=Tsit5(), kwargs...)

Return the effective temporal mode required for input cavity `i` when several source
cavities are cascaded before the physical system. The other source cavities modify the field
propagated through the input chain, so the nominal mode `u_fcts[i]` must be corrected before
its coupling is constructed.

The input modes must be ordered starting with the first input cavity before the system, and
`i > 1`. `u_fcts` may contain callables or sampled mode vectors. When `gu_fcts` is omitted,
[`coupling_input`](@ref) is used to construct the couplings of the nominal modes. Sampled
coupling data may also be supplied explicitly.

Returns a callable `uᵢ_eff(t)`. `alg` and all remaining keyword arguments are forwarded to
the auxiliary OrdinaryDiffEq solve.

See also [`effective_output_mode`](@ref) and [Pulse modes](@ref).
"""
function effective_input_mode(u_fcts, gu_fcts, T, i; alg = Tsit5(), kwargs...)
    @assert i > 1
    n = i - 1
    # Capture as a tuple for concrete closure types
    gu = ntuple(k -> gu_fcts[k], n)
    u_i = u_fcts[i]
    gu_i = gu_fcts[i]
    function multiple_inputs_α!(dα, α, p, t)
        gu_buf = p
        @inbounds for k = 1:n
            gu_buf[k] = gu[k](t)
        end
        ui_t = u_i(t)
        gui_t = gu_i(t)
        @inbounds for j = 1:n
            coupling_sum = zero(ComplexF64)
            for k = 1:(j-1)
                coupling_sum += gui_t' * α[k]
            end
            dα[j] = -gu_buf[j] * (ui_t - coupling_sum) + 0.5 * abs2(gu_buf[j]) * α[j]
        end
    end
    u0 = zeros(ComplexF64, n)
    p = Vector{ComplexF64}(undef, n)
    tspan = (T[1], T[end])
    prob = ODEProblem(multiple_inputs_α!, u0, tspan, p)
    sol_α = solve(prob, alg; kwargs...)
    function u_i_eff(t)
        α_t = sol_α(t)
        result = u_i(t)
        @inbounds for k = 1:n
            result -= gu[k](t)' * α_t[k]
        end
        return result
    end
    return u_i_eff
end

effective_input_mode(u_fcts, T, i; kwargs...) =
    effective_input_mode(u_fcts, [coupling_input(u_, T) for u_ in u_fcts], T, i; kwargs...)

function effective_input_mode(
    u_data::AbstractVector{<:AbstractVector},
    gu_data::AbstractVector{<:AbstractVector},
    T::AbstractVector,
    i;
    alg = Tsit5(),
    kwargs...,
)
    u_fcts = [_interpolate_mode(u_, T) for u_ in u_data]
    gu_fcts = [_interpolate_mode(gu_, T) for gu_ in gu_data]
    return effective_input_mode(u_fcts, gu_fcts, T, i; alg, kwargs...)
end

effective_input_mode(
    u_data::AbstractVector{<:AbstractVector},
    T::AbstractVector,
    i;
    alg = Tsit5(),
    kwargs...,
) = effective_input_mode(
    u_data,
    [coupling_input(u_, T).(T) for u_ in u_data],
    T,
    i;
    alg,
    kwargs...,
)

# ──────────────────────────────────────────────
# coupling_delay_out / coupling_delay_in
# ──────────────────────────────────────────────

function _delay_coupling_from_modes(num_mode::Vector, u::Vector, v::Vector, T::Vector)
    nt = length(T)
    mode_sq = Vector{Float64}(undef, nt)
    map!(abs2, mode_sq, u)
    ∫u2 = cumul_integrate(T, mode_sq)
    map!(abs2, mode_sq, v)
    ∫v2 = cumul_integrate(T, mode_sq)
    coupling = zeros(ComplexF64, nt)
    @inbounds for i = 1:nt
        denom = abs(∫v2[i] - ∫u2[i])
        if sqrt(abs(denom)) > _tol_div
            coupling[i] = num_mode[i]' / sqrt(denom + _ϵ)
        end
    end
    return PulseCoupling(LinearInterpolation(coupling, T; extrapolation = _extrapolate))
end

"""
    coupling_delay_out(u, v, T)

Construct the emission coupling of a virtual delay cavity that simultaneously emits the
mode `u` and absorbs the mode `v`. For

```math
D(t)=\int_{T[1]}^t \left(|v(t')|^2-|u(t')|^2\right)dt',
```

the ideal coupling is

```math
\tilde g_{\mathrm{out}}(t)=\frac{u^*(t)}{\sqrt{D(t)}}.
```

`u` and `v` may be sampled vectors, functions, or linear interpolations. The physical delay
construction requires `D(t) ≥ 0`; the sampled implementation evaluates the denominator as
`√(abs(D) + ϵ)` and returns zero when it is below the division tolerance. The result is a
[`PulseCoupling`](@ref).

Use together with [`coupling_delay_in`](@ref). See [Quantum pulses and input-output theory](@ref)
for the delay-cavity interpretation.
"""
coupling_delay_out(u::Vector, v::Vector, T::Vector) = _delay_coupling_from_modes(u, u, v, T)
coupling_delay_out(u::Function, v::Function, T::Vector) =
    coupling_delay_out(u.(T), v.(T), T)
coupling_delay_out(u::LinearInterpolation, v::LinearInterpolation, T::Vector) =
    coupling_delay_out(u.(T), v.(T), T)

"""
    coupling_delay_in(u, v, T)

Construct the absorption coupling of a virtual delay cavity that simultaneously emits the
mode `u` and absorbs the mode `v`. With the same stored norm `D(t)` used by
[`coupling_delay_out`](@ref), the ideal coupling is

```math
\tilde g_{\mathrm{in}}(t)=-\frac{v^*(t)}{\sqrt{D(t)}}.
```

`u` and `v` may be sampled vectors, functions, or linear interpolations. The same physical
condition `D(t) ≥ 0` and sampled denominator handling as in [`coupling_delay_out`](@ref)
applies. The result is a [`PulseCoupling`](@ref).
"""
coupling_delay_in(u::Vector, v::Vector, T::Vector) = _delay_coupling_from_modes(-v, u, v, T)
coupling_delay_in(u::Function, v::Function, T::Vector) = coupling_delay_in(u.(T), v.(T), T)
coupling_delay_in(u::LinearInterpolation, v::LinearInterpolation, T::Vector) =
    coupling_delay_in(u.(T), v.(T), T)
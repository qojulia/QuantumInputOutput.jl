```@meta
CurrentModule = QuantumInputOutput
CollapsedDocStrings = true
```

# Pulse modes

Traveling wave packets live in a continuum, whereas the open-system solvers used by QuantumInputOutput evolve discrete quantum degrees of freedom. The virtual-cavity construction bridges the two: a chosen temporal mode is represented by an auxiliary oscillator whose time-dependent coupling emits into, or absorbs from, that mode [Kiilerich2019, Kiilerich2020](@cite).

The theory behind the construction is given in [Quantum pulses and input-output theory](@ref). This page is about assembling pulse modes in a model.

## Single input and output modes

A temporal mode can be supplied as sampled data or as a function evaluated on a time grid. Gaussian wave packets also have an analytic representation.

```@docs
Gaussian
PulseCoupling
```

Place a virtual source cavity upstream of the physical system and a virtual receiver cavity downstream.

```@docs
coupling_input
coupling_output
```

A typical symbolic model contains scalar coupling parameters whose time dependence is supplied only when the model is translated:

```julia
G_u = SLH(1, g_u * a_u, 0)
G_s = SLH(1, L_s, H_s)
G_v = SLH(1, g_v * a_v, 0)
G = G_u ▷ G_s ▷ G_v

T = collect(range(0, 12; length = 4001))
gu = coupling_input(Gaussian(4.0, 1.0))
gv = coupling_output(v, T)
```

## Several selected modes

Several explicit virtual cavities alter the fields seen by one another. Their ordering therefore carries physical information: input modes follow the source-cavity cascade before the system, while output modes follow the receiver-cavity cascade after it.

```@docs
effective_input_mode
effective_output_mode
```

The corrected mode is used to construct the corresponding single-mode coupling:

```julia
u2_eff = effective_input_mode(u_modes, T, 2)
gu2 = coupling_input(u2_eff, T)

v2_eff = effective_output_mode(v_modes, T, 2)
gv2 = coupling_output(v2_eff, T)
```

## Propagation delays

A short delay can require a virtual cavity to absorb one envelope while emitting another at the same time.

```@docs
coupling_delay_in
coupling_delay_out
```

The delay construction follows [Christiansen2026](@cite) and is derived in [Quantum pulses and input-output theory](@ref). See the [pulse-delay example](../examples/08-1_pulse-delay__simple.md) for a complete network.

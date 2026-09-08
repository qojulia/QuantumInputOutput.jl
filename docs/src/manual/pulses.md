```@meta
CurrentModule = QuantumInputOutput
CollapsedDocStrings = true
```

# Pulse modes

Traveling wave packets live in a continuum, whereas the open-system solvers used by QuantumInputOutput evolve discrete quantum degrees of freedom. The virtual-cavity construction bridges the two: a chosen temporal mode is represented by an auxiliary oscillator whose time-dependent coupling emits into, or absorbs from, that mode [Kiilerich2019, Kiilerich2020](@cite).

The theory behind the construction is given in [Quantum pulses and input-output theory](@ref). This page is about assembling pulse modes in a model.

## Single input and output modes

A temporal mode can be supplied as sampled data or as a function evaluated on a time grid. For Gaussian wave packets, `Gaussian` provides an analytic route that avoids numerical integration of the mode envelope.

```@docs
Gaussian
```

The coupling returned by the pulse constructors is a callable wrapper designed to be used directly as a time-dependent scalar parameter during numerical translation.

```@docs
PulseCoupling
```

Use an input coupling for a virtual cavity placed upstream of the physical system and an output coupling for a virtual cavity placed downstream. The two constructions are directional: the first releases an initially occupied auxiliary mode into the field, while the second captures a prescribed outgoing temporal mode.

```@docs
coupling_input
coupling_output
```

A typical symbolic model therefore contains scalar parameters such as `g_u` and `g_v`, while their numerical values are supplied only at translation time:

```julia
G_u = SLH(1, g_u * a_u, 0)
G_s = SLH(1, L_s, H_s)
G_v = SLH(1, g_v * a_v, 0)
G = G_u ▷ G_s ▷ G_v

T = range(0, 12; length = 4001)
gu = coupling_input(Gaussian(4.0, 1.0))
gv = coupling_output(v, collect(T))
```

## Several selected modes

Virtual cavities themselves scatter the fields seen by later cavities. Consequently, when several input or output modes are represented explicitly, the nominal temporal modes cannot in general be coupled independently. The effective-mode helpers account for that distortion before the corresponding coupling is constructed.

The ordering of the mode collection is physical. For inputs it follows the sequence of virtual input cavities before the system; for outputs it follows the sequence of virtual output cavities after the system.

```@docs
effective_input_mode
effective_output_mode
```

The usual pattern is to determine the effective mode first and then construct its coupling:

```julia
u2_eff = effective_input_mode(u_modes, T, 2)
gu2 = coupling_input(u2_eff, T)

v2_eff = effective_output_mode(v_modes, T, 2)
gv2 = coupling_output(v2_eff, T)
```

## Propagation delays

A finite delay can also be represented by a virtual cavity. During a short delay the auxiliary mode may need to absorb one envelope while simultaneously emitting another, so the ordinary source or receiver coupling is not sufficient. `coupling_delay_in` and `coupling_delay_out` construct the two couplings for that case.

```@docs
coupling_delay_in
coupling_delay_out
```

The delay construction is discussed in [Quantum pulses and input-output theory](@ref) and follows the finite-propagation treatment of [Christiansen2026](@cite). See [Pulse delay](@ref) in the examples for a complete network.

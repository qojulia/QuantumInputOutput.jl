```@meta
CurrentModule = QuantumInputOutput
CollapsedDocStrings = true
```

# SLH networks

QuantumInputOutput represents a component by its scattering, coupling, and Hamiltonian data and keeps that representation compositional. The same network algebra can be used with symbolic operators from [SecondQuantizedAlgebra.jl](@extref SecondQuantizedAlgebra :doc:`index`) or with numerical QuantumOptics.jl operators; the distinction matters only when the finished model is handed to a numerical backend.

## Components

An ``N``-port component has an ``N\times N`` scattering matrix, one coupling operator per port, and one Hamiltonian. Port order is part of the model: composition connects channels by their positions in these arrays.

```@docs
SLH
```

The three pieces remain directly accessible after composition. This is usually the point at which a symbolic network is inspected, transformed, or translated.

```@docs
scattering
jump_operator
hamiltonian
```

`lindblad` is retained only as a deprecated name for `jump_operator`.

```@docs
lindblad
```

## Composition

Network topology is written algebraically. Cascade connects equally sized channel sets in series; concatenation places independent components beside one another. The Unicode operators are convenient in equations, while the named forms are useful in ordinary code and accept the same variadic composition.

```@docs
▷
cascade
```

```@docs
⊞
concatenate
```

Internal coherent loops are eliminated at the SLH level. When several connections are supplied together, their port labels refer to the original unreduced component, so the network can be specified without manually renumbering ports after each reduction.

```@docs
feedback
```

A useful modeling discipline is therefore to build the physical network first and inspect only the final ``(S,L,H)`` triple. The intermediate algebra is an implementation of the network topology, not a sequence of master equations that has to be derived by hand.

## From a network to dynamics

For a symbolic model, numerical translation is best delayed until composition and feedback reduction are complete. `to_numeric` lowers the Hamiltonian and coupling operators onto a QuantumOptics.jl basis while preserving time-dependent parameters such as pulse couplings.

```@docs
to_numeric
```

QuantumInputOutput does not prescribe the subsequent solver. Full Hilbert-space dynamics can be evolved with QuantumOptics.jl, while the symbolic Hamiltonian and coupling operators can also be passed to moment-based tools such as QuantumCumulants.jl. The complete full-Hilbert-space workflow is shown in the [Tutorial](@ref).

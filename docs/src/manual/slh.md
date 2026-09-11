```@meta
CurrentModule = QuantumInputOutput
CollapsedDocStrings = true
```

# SLH networks

QuantumInputOutput keeps the physical network and the numerical solver as separate concerns. Components are represented in SLH form and composed at that level; only the resulting model is translated to the representation used for dynamics. The same network algebra works with symbolic operators from [SecondQuantizedAlgebra.jl](@extref SecondQuantizedAlgebra :doc:`index`) and with numerical QuantumOptics.jl operators.

## Components

Port order is part of the network topology: channel positions determine which fields are connected by later composition operations.

```@docs
SLH
```

The components of a composed model remain accessible independently.

```@docs
scattering
jump_operator
hamiltonian
lindblad
```

## Composition

Cascade describes propagation through successive components, while concatenation places independent channels beside one another. These operations keep the model in SLH form, so network structure can be changed without re-deriving a master equation.

```@docs
▷
cascade
⊞
concatenate
```

Feedback closes internal coherent connections before numerical translation.

```@docs
feedback
```

A useful modeling discipline is therefore to construct the physical network completely and inspect only the final ``(S,L,H)`` triple.

## From a network to dynamics

Symbolic construction is useful when the composed model must still be inspected, simplified, transformed, or reused by different backends. If none of those operations is needed, numerical QuantumOptics.jl operators and time-dependent callables can instead be placed directly in `SLH` objects, avoiding a separate symbolic-lowering step.

For symbolic models, translation is normally the boundary between model construction and numerical evolution.

```@docs
to_numeric
```

The symbolic and numerical Hilbert spaces do not have to be identical. Operator substitutions forwarded through `to_numeric(...; operators=...)` can map selected symbolic operators onto a reduced numerical basis—for example when an auxiliary output cavity was useful in the symbolic derivation but is intentionally omitted from a particular solve.

QuantumInputOutput does not prescribe the subsequent solver. Full Hilbert-space dynamics can be evolved with QuantumOptics.jl, while the symbolic Hamiltonian and jump operators can also be passed to moment-based tools such as QuantumCumulants.jl. The complete full-Hilbert-space workflow is shown in the [Tutorial](@ref).

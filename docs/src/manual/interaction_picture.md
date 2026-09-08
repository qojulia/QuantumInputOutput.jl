```@meta
CurrentModule = QuantumInputOutput
CollapsedDocStrings = true
```

# Interaction picture

Virtual input and output cavities can exchange a large number of excitations even when the physically interesting subsystem occupies only a small part of its Hilbert space. For high-photon-number pulse problems it can therefore be advantageous to transform away the deterministic transfer among the selected pulse modes and simulate only the residual dynamics [Christiansen2023](@cite).

The derivation is given in [Interaction-picture formulation](@ref). At the API level the construction has two steps: build the mode-transfer generator, then integrate the corresponding change of basis.

## Coupling matrix

Supply the virtual-mode couplings in their physical cascade order.

```@docs
coupling_matrix
```

## Mode evolution

The resulting transformation is evaluated when substituting the virtual-mode operators into the interaction-picture model.

```@docs
solve_mode_evolution
```

For the symmetric source/receiver case, the same transformation is available in closed form.

```@docs
solve_mode_evolution_symmetric
```

A minimal setup is

```julia
A = coupling_matrix((gu, gv))
M = solve_mode_evolution(A, T)

# M(t) maps the original virtual-mode operators into the interaction picture.
Mt = M(T[end] / 2)
```

The package deliberately returns the mode transformation rather than rewriting a complete SLH model automatically. This keeps the representation change explicit: construct the network first, identify the pure virtual-mode transfer to remove, transform the affected operators, and then translate the resulting symbolic model to the desired numerical backend.

See the [interaction-picture scattering example](../examples/06-1_interaction-picture__PRA2023_107-013706_fig2.md) for the complete high-photon-number workflow.

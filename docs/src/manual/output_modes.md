```@meta
CurrentModule = QuantumInputOutput
CollapsedDocStrings = true
```

# Output modes

The outgoing field of a nonlinear quantum system is generally not known in advance and need not occupy the same temporal mode as the input. QuantumInputOutput therefore separates two tasks:

1. evolve the system while leaving the relevant output continuum unprojected;
2. identify the populated temporal modes from its first-order correlation kernel and, when needed, add selected modes back as explicit virtual output cavities.

This closes the loop from a discrete input state to a traveling output field and back to a discrete output-mode state.

## Correlation kernel

The first-order correlation kernel is the bridge between the continuous output field and a discrete temporal-mode basis.

```@docs
correlation_matrix
```

## Extracting temporal modes

Diagonalizing the Hermitian sampled kernel gives the natural temporal-mode basis. For a uniform grid with spacing ``\Delta t``, the discrete eigenvectors are converted to continuum-normalized modes by dividing by ``\sqrt{\Delta t}``, while the photon occupation associated with an eigenvalue ``\lambda_i`` is ``n_i \approx \lambda_i\Delta t``.

```julia
G1 = correlation_matrix(T, ρt, H, J, Ls)
F = eigen(G1)

Δt = T[2] - T[1]
occupations = real.(F.values) .* Δt
v1 = F.vectors[:, end] ./ sqrt(Δt)
```

The eigenvectors are defined only up to an overall phase. If a later calculation depends on a particular phase convention, fix that phase explicitly before constructing the virtual-cavity coupling.

## Capturing a selected output mode

Once a relevant mode has been identified, use [`coupling_output`](@ref) and repeat the evolution with a downstream virtual cavity. The auxiliary cavity then stores the quantum state associated with that temporal mode, so ordinary discrete-mode observables and state diagnostics can be applied.

```julia
gv = coupling_output(v1, T)
```

For several simultaneously retained output modes, first account for the disturbance introduced by preceding receiver cavities with [`effective_output_mode`](@ref). The full single-photon workflow is carried out in the [Tutorial](@ref), while the continuum-mode interpretation is developed in [Quantum pulses and input-output theory](@ref).

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

The continuous modes solve an integral eigenvalue problem, so a sampled kernel should be paired with the quadrature rule used to approximate that integral. Given positive quadrature weights ``w_j`` and ``W=\operatorname{diag}(w_j)``, diagonalize the Hermitian matrix ``W^{1/2} K W^{1/2}``. Its eigenvalues are the mode occupations directly, while the corresponding continuum-normalized mode samples are obtained by applying ``W^{-1/2}`` to its eigenvectors.

```julia
K = correlation_matrix(T, ρt, H, J, Ls)

sqrtw = sqrt.(w)
Kweighted = Hermitian(sqrtw .* Matrix(K) .* transpose(sqrtw))
F = eigen(Kweighted)

occupations = real.(F.values)
v1 = F.vectors[:, end] ./ sqrtw
```

For the rectangular rule on a uniform grid, ``w_j=\Delta t`` for every sample. Then ``W`` is proportional to the identity, so the shortcut

```julia
F = eigen(K)
Δt = T[2] - T[1]
occupations = real.(F.values) .* Δt
v1 = F.vectors[:, end] ./ sqrt(Δt)
```

is equivalent. Other rules—including trapezoidal quadrature on a uniform grid—should use their actual weights in the weighted construction above.

The eigenvectors are defined only up to an overall phase. If a later calculation depends on a particular phase convention, fix that phase explicitly before constructing the virtual-cavity coupling.

## Capturing a selected output mode

Once a relevant mode has been identified, use [`coupling_output`](@ref) and repeat the evolution with a downstream virtual cavity. The auxiliary cavity then stores the quantum state associated with that temporal mode, so ordinary discrete-mode observables and state diagnostics can be applied.

```julia
gv = coupling_output(v1, T)
```

For several simultaneously retained output modes, first account for the disturbance introduced by preceding receiver cavities with [`effective_output_mode`](@ref). The full single-photon workflow is carried out in the [Tutorial](@ref), while the continuum-mode interpretation is developed in [Quantum pulses and input-output theory](@ref).

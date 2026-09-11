```@raw html
---
layout: home

hero:
  name: QuantumInputOutput.jl
  text: Composable input-output models for propagating quantum pulses
  tagline: Write the quantum network, not the master equation.
  actions:
    - theme: brand
      text: Tutorial
      link: tutorial.md
    - theme: alt
      text: Theory
      link: theory.md
    - theme: alt
      text: View on GitHub
      link: https://github.com/qojulia/QuantumInputOutput.jl
  image:
    src: /assets/logo.svg
    alt: QuantumInputOutput.jl Gaussian-pulse logo

features:
  - icon: ⇢
    title: Compose quantum networks
    details: Build optical networks from SLH components using cascade, concatenation, phase shifts, and feedback reduction.
  - icon: 〰️
    title: Embed traveling pulses
    details: Represent selected temporal input and output modes as virtual cavities inside an ordinary open quantum system.
  - icon: ↓
    title: Lower to numerical solvers
    details: Translate symbolic SLH models into numerical representations for full Hilbert-space or cumulant-based dynamics.
  - icon: ∫
    title: Recover output modes
    details: Compute field correlations, identify populated temporal modes, and promote selected outgoing modes back into explicit quantum subsystems.
---
```

```@meta
CurrentModule = QuantumInputOutput
```

`QuantumInputOutput.jl` is a Julia framework for constructing and simulating quantum input-output networks with propagating pulses. The package separates **physical model construction** from **numerical evolution**: users describe the network and the traveling modes, while QuantumInputOutput derives the corresponding symbolic Hamiltonian, coupling operators, and scattering structure needed by standard open-system solvers.

A typical workflow is

```text
physical network + propagating pulse
                ↓
        compositional SLH model
                ↓
        symbolic (S, L, H)
                ↓
          numerical backend
                ↓
 dynamics + outgoing temporal modes
```

## Get Started

Install the package with Julia's package manager:

```julia-repl
pkg> add QuantumInputOutput
```

Then start with the [Tutorial](tutorial.md), which walks through the complete pulse-scattering workflow from model construction to output-mode extraction.

## Core Workflow

### Compose the network

Represent optical components as SLH elements and combine them algebraically using cascade, concatenation, phase elements, and feedback. This keeps the model close to the physical network instead of requiring a manual derivation of the final master equation.

### Embed propagating pulses

Selected temporal input and output modes are represented by virtual cavities with time-dependent couplings. Traveling quantum states therefore become ordinary quantum degrees of freedom that can be evolved with standard open-system methods.

### Choose the numerical representation

The symbolic model can be translated to numerical operators for full Hilbert-space calculations with QuantumOptics.jl, or used with QuantumCumulants.jl for reduced moment dynamics. QuantumInputOutput owns the network and pulse representation rather than a single numerical solver.

### Analyze the outgoing field

Two-time field correlations provide a temporal-mode decomposition of the outgoing radiation. Dominant output modes can then be treated as explicit quantum subsystems for state and observable calculations.

See the [Theory](theory.md) for the underlying input-output and virtual-cavity formalism, [Implementation](implementation.md) for the symbolic-to-numeric pipeline, [Reference](api.md) for the API, and [Relevant Literature](literature.md) for the foundational papers.

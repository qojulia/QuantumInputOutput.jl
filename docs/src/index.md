```@raw html
---
layout: home

hero:
  name: QuantumInputOutput.jl
  text: Composable input-output models for propagating quantum pulses
  tagline: Symbolic SLH network construction and numerical pulse dynamics in Julia.
  actions:
    - theme: brand
      text: Tutorial
      link: tutorial.md
    - theme: alt
      text: Theory
      link: theory.md
    - theme: alt
      text: Implementation
      link: implementation.md
    - theme: alt
      text: View on GitHub
      link: https://github.com/qojulia/QuantumInputOutput.jl
  image:
    src: /assets/logo.svg
    alt: QuantumInputOutput.jl Gaussian-pulse logo

features:
  - icon: ⇢
    title: Compose quantum networks
    details: Build quantum networks from SLH components using cascade, concatenation, and feedback reduction.
  - icon: 〰️
    title: Embed traveling pulses
    details: Represent selected temporal input and output modes as virtual cavities inside an ordinary open quantum system.
  - icon: ℏ
    title: Translate to numerics
    details: Translate symbolic SLH models into numerical representations for full Hilbert-space or cumulant-based dynamics.
  - icon: ∫
    title: Recover output modes
    details: Compute field correlations, identify populated temporal modes, and promote selected outgoing modes back into explicit quantum subsystems.
  - icon: Δt
    title: Extend pulse models
    details: Use interaction-picture transformations for efficient simulations and virtual delay cavities to model pulse propagation delays.
---
```

```@meta
CurrentModule = QuantumInputOutput
```

`QuantumInputOutput.jl` is a Julia framework for constructing and simulating quantum input-output networks with propagating pulses. The package separates **physical model construction** from **numerical evolution**: users describe the network and the traveling modes, while QuantumInputOutput derives the corresponding symbolic Hamiltonian, coupling operators, and scattering structure needed by standard open-system solvers.

## Core Workflow

### 1. Compose the network

Represent quantum systems and their input-output connections as SLH elements, then combine them algebraically using cascade, concatenation, and feedback reduction. This keeps the model close to the physical network instead of requiring a manual derivation of the final master equation.

### 2. Embed propagating pulses

Selected temporal input and output modes are represented by virtual cavities with time-dependent couplings. Traveling quantum states therefore become ordinary quantum degrees of freedom that can be evolved with standard open-system methods.

### 3. Choose the numerical representation

The symbolic model can be translated to numerical operators for full Hilbert-space calculations with QuantumOptics.jl, or used with QuantumCumulants.jl for reduced moment dynamics. QuantumInputOutput owns the network and pulse representation rather than a single numerical solver.

### 4. Analyze the outgoing field

Two-time field correlations provide a temporal-mode decomposition of the outgoing radiation. Dominant output modes can then be treated as explicit quantum subsystems for state and observable calculations.

See the [Theory](theory.md) for the underlying input-output and virtual-cavity formalism, [Implementation](implementation.md) for the symbolic-to-numeric pipeline, [API](api.md) for the public interface, and [Relevant Literature](literature.md) for the foundational papers.

## Get Started

Install the package with Julia's package manager:

```julia-repl
pkg> add QuantumInputOutput
```

Then start with the [Tutorial](tutorial.md), which walks through the complete pulse-scattering workflow from model construction to output-mode extraction.

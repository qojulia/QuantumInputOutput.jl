###########################
### interaction picture ###
###########################

using StaticArrays: SMatrix

_as_time_function(x::Number) = _ -> x
_as_time_function(x) = x  # anything callable passes through

"""
    coupling_matrix(gs::Tuple)
    coupling_matrix(g1, g2, gs...)

Construct the anti-Hermitian coefficient matrix `A(t)` that generates the passive
interaction-picture mixing of a cascade of virtual modes. The couplings must be supplied
in physical cascade order and may be constants or callables `gᵢ(t)`.

```math
A_{ij}(t) = \frac{1}{2}
\begin{cases}
0, & i=j, \\
g_i(t)g_j^*(t), & i<j, \\
-g_j^*(t)g_i(t), & i>j.
\end{cases}
```

Returns a closure `t -> A(t)` whose value is an `N × N` static `ComplexF64` matrix.
Because `A(t)' == -A(t)`, the solution of `dM/dt = A(t)M(t)` is unitary.

See also [`solve_mode_evolution`](@ref) and [Interaction-picture formulation](@ref).
"""
function coupling_matrix(gs::Tuple{Vararg{Any,N}}) where {N}
    gfs = map(_as_time_function, gs)
    A(t) = _coupling_matrix(gfs, t)
    return A
end

# Unrolled at compile time so that tuples of mixed element type (a constant next to
# an interpolant, say) stay allocation-free: indexing such a tuple in a loop is only
# type stable if constant propagation reaches the index, and it does not.
@generated function _coupling_matrix(gfs::Tuple{Vararg{Any,N}}, t) where {N}
    g = [Symbol(:g_, i) for i = 1:N]
    calls = [:($(g[i]) = ComplexF64(gfs[$i](t))) for i = 1:N]
    entries = Expr[]
    for j = 1:N, i = 1:N  # column-major
        e = if i == j
            :(zero(ComplexF64))
        elseif i < j
            :(0.5 * $(g[i]) * conj($(g[j])))
        else
            :(-0.5 * conj($(g[j])) * $(g[i]))
        end
        push!(entries, e)
    end
    quote
        $(calls...)
        SMatrix{$N,$N,ComplexF64,$(N * N)}(($(entries...),))
    end
end

coupling_matrix(g1, g2, gs...) = coupling_matrix((g1, g2, gs...))

"""
    solve_mode_evolution(A::Function, T; alg=Tsit5(), kwargs...)

Solve the interaction-picture matrix equation

```math
\dot M(t)=A(t)M(t), \qquad M(T[1])=I,
```

and save the solution on the grid `T`. The matrix size is inferred from `A(T[1])`.
`A` is normally constructed with [`coupling_matrix`](@ref).

Returns the OrdinaryDiffEq solution directly, so `sol(t)` evaluates the mode-transformation
matrix at arbitrary times supported by the solver interpolation. `alg` selects the ODE
algorithm and all remaining keyword arguments are forwarded to `solve`.
"""
function solve_mode_evolution(A::Function, T; alg = Tsit5(), kwargs...)
    T0 = T[1]
    Tend = T[end]
    n = size(A(T0), 1)
    M0 = Matrix{ComplexF64}(I, n, n)
    function f_M!(du, u, p, t)
        mul!(du, A(t), u)
    end
    prob = ODEProblem(f_M!, M0, (T0, Tend))
    sol = solve(prob, alg; saveat = T, kwargs...)
    return sol
end

"""
    solve_mode_evolution_symmetric(u, T)

Return the analytic two-mode interaction-picture rotation for a source and receiver built
from the same normalized temporal envelope `u`. `u` may be a callable or values sampled on
`T`.

With

```math
F(t)=\int_{T[1]}^t |u(t')|^2\,dt', \qquad \sin^2\theta(t)=F(t),
```

the returned callable evaluates

```math
M(t)=\begin{bmatrix}
\cos\theta(t) & -\sin\theta(t) \\
\sin\theta(t) & \cos\theta(t)
\end{bmatrix}.
```

The accumulated norm is clamped to `[0, 1]` before evaluating the angle, which suppresses
small numerical excursions outside the physical interval. Matrix elements are linearly
interpolated over `T` and extrapolated using the boundary extension.

See also [`solve_mode_evolution`](@ref).
"""
function solve_mode_evolution_symmetric(u, T)
    u_vals = u isa Function ? u.(T) : u
    sin2θ = cumul_integrate(T, abs2.(u_vals))
    sin2θ = clamp.(real.(sin2θ), 0.0, 1.0)
    θ = asin.(sqrt.(sin2θ))
    cθ = cos.(θ)
    sθ = sin.(θ)

    M11 = LinearInterpolation(cθ, T; extrapolation = ExtrapolationType.Extension)
    M12 = LinearInterpolation(-sθ, T; extrapolation = ExtrapolationType.Extension)
    M21 = LinearInterpolation(sθ, T; extrapolation = ExtrapolationType.Extension)
    M22 = LinearInterpolation(cθ, T; extrapolation = ExtrapolationType.Extension)

    return t -> SMatrix{2,2}(M11(t), M21(t), M12(t), M22(t))  # column-major
end

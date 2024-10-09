# SplitODEProblem

The `SplitODEProblem` type from package [DifferentialEquations.jl](https://diffeq.sciml.ai/stable/types/split_ode_types/) offers a interface for problems similar to the ones we are trying to solve.

Consider the Henon-Heiles system we can use one of the solvers dedicated to split problems, see [Split ODE Solvers](https://diffeq.sciml.ai/stable/solvers/split_ode_solve/#split_ode_solve):

```@example 4
using Plots
```

```@example 4
using OrdinaryDiffEq

epsilon = 0.002
A = [ 0 0 1 0 ;
      0 0 0 0 ;
     -1 0 0 0 ;
      0 0 0 0 ]

f1 = DiffEqArrayOperator( A ./ epsilon)

function f2(du, u, p, t)
    du[1] = 0
    du[2] = u[4]
    du[3] = 2*u[1]*u[2]
    du[4] = -u[2] - u[1]^2 + u[2]^2 
end

tspan = (0.0, 0.1)

u0 = [0.55, 0.12, 0.03, 0.89]

prob1 = SplitODEProblem(f1, f2, u0, tspan);
sol1 = solve(prob1, ETDRK4(), dt=0.001);
nothing # hide
```
With our method we need to give the value of `epsilon`. In some case, you can get this value from `f1`.
You can pass a `DiffEqArrayOperator` as argument to the problem but the method is not always valid 
so we define a new type called `LinearHOODEOperator`:

```@example 4
using HOODESolver

f1 = LinearHOODEOperator(epsilon, A)
prob2 = SplitODEProblem(f1, f2, u0, tspan)
sol2 = solve(prob2, HOODEAB(), dt=0.01)

plot(sol1, vars=[3], label="EDTRK4")
plot!(sol2, vars=[3], label="HOODEAB")
plot!(sol2.t, getindex.(sol2.u, 3), m=:o)
```

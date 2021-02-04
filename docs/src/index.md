# HOODESolver.jl

The objective of this Julia package is to valorize the recent developments carried out within [INRIA team MINGuS](https://team.inria.fr/mingus/) on Uniformly Accurate numerical methods (UA) for highly oscillating problems. We propose to solve the following equation 

$$\frac{d u(t)}{dt} = \frac{1}{\varepsilon} A u(t) + f(t, u(t)), \qquad u(t=t_{start})=u_{in}, \qquad \varepsilon\in ]0, 1], \qquad (1)$$

with 
-  $u : t\in [t_{start}, t_{end}] \mapsto u(t)\in \mathbb{R}^n, \quad t_{start}, t_{end}\in \mathbb{R}$, 
-  $u_{in}\in \mathbb{R}^n$, 
-  $A\in {\mathcal{M}}_{n,n}(\mathbb{R})$ is such that $\tau \mapsto \exp(\tau A)$ is 2 \pi -periodic,  
-  $f : (t, u) \in  \mathbb{R}\times \mathbb{R}^n \mapsto \mathbb{R}^n$.

The purpose here is to write an explanatory documentation of the *Julia* package containing the two-scale method (see [chartier2020](@cite), [chartier2015](@cite) and [crouseilles2013](@cite). This package is inspired by the Differential Equations package [SciML](https://diffeq.sciml.ai/dev/index.html).


## References

```@bibliography
```

## Index

```@index
```

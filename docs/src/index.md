# HOODESolver.jl

The objective of this Julia package is to valorize the recent developments carried out within MINGuS on Uniformly Accurate numerical methods (UA) for highly oscillating problems. We propose to solve the following equation 

$$\frac{d u(t)}{dt} = \frac{1}{\varepsilon} A u(t) + f(u(t)), \;\;\; u(t=t_{0})=u_{0}, \;\; \varepsilon\in ]0, 1], \;\;\;\;\;\;\;\;\;\;\;\; (1)$$

with 
-  $u : t\in [t_{0}, t_{fin}] \mapsto u(t)\in \mathbb{R}^n, \;\; t_{0}, t_{\text{end}}\in \mathbb{R}$, 
-  $u_{0}\in \mathbb{R}^n$, 
-  $A\in {\mathcal{M}}_{n,n}(\mathbb{R})$ is such that $\tau \mapsto \exp(\tau A)$ is periodic,  
-  $f : u\in  \mathbb{R}^n \mapsto \mathbb{R}^n$.

The purpose here is to write an explanatory documentation of the *Julia* package containing the two-scale method (see [chartier2020](@cite), [chartier2015](@cite) and [crouseilles2013](@cite). This package is inspired by the Differential Equations package [SciML](https://diffeq.sciml.ai/dev/index.html).


## References

```@bibliography
```

## Index

```@index
```

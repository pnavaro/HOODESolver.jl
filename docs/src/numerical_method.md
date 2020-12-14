
# Numerical method

## Two-scale formulation 

First, rewrite equation (1) using the variable change $w(t)=\exp(-(t-t_{0})A/\varepsilon) u(t)$ to obtain

$$\frac{d w(t)}{dt} = F\Big( \frac{t-t_{0}}{\varepsilon}, w(t) \Big), \;\;\; w(t_{0})=u_{0}, \;\; \varepsilon\in ]0, 1],$$

where the function $F$ is expressed from the data of the original problem (1)

$$F\Big( \frac{s}{\varepsilon}, w \Big) = \exp(-sA/\varepsilon) \; f( \exp(sA/\varepsilon), \; w).$$

We then introduce the function $U(t, \tau), \tau\in [0, 2 \pi]$ such that $U(t, \tau=(t-t_{0})/\varepsilon) = w(t)$. The two-scale function is then the solution of the following equation

$$\frac{\partial U}{\partial t} + \frac{1}{\varepsilon} \frac{\partial U}{\partial \tau} =  F( \tau, U), \;\;\; U(t=t_{0}, \tau)=\Phi(\tau), \;\; \varepsilon\in ]0, 1], \;\;\;\;\;\;\;\;\;\; (2)$$

where $\Phi$ is a function checking $\Phi(\tau=0)=u_{0}$ chosen so that the $U$ solution of (2) is regular (see [chartier2015](@cite) and [chartier2020](@cite).

## Discretization  

The numerical method is based on a discretization of equation (2). In the direction $\tau$, a spectral method is used, while for the time $t$, an exponential Adams-Bashforth method allows to build a high order method (see [chartier2020](@cite)). The initialization is based on a "butterfly" technique (going back and forth around the initial time).

### Initialization 

Let r be the order of the method $AB_r$.\
Let $\Delta t$ the time step, for $i \in \{r, -(r-1), \ldots, r-1, r\}$,
we note $u_i = u(t_0+i \Delta t)$.\
Let $r'$ be the orders of the intermediate AB methods we will use.\
If $u_{k}$ is known with a precision of ${\mathcal O}(\Delta t^{r'+1})$, and for $r' \geq 2, u_{k-1}, \ldots, u_{k-r'+1}$ are known with a precision of ${\mathcal O}(\Delta t^{r'})$ then we can calculate $u_{k+1}$ with a precision of ${\mathcal O}(\Delta t^{r'+1})$ with the method $AB_{r'}$.\
Similarly, if $u_{k}$ is known with a precision of ${\mathcal O}(\Delta t^{r'+1})$, and for $r' \geq 2, u_{k+1}, \ldots, u_{k+r'-1}$ are known with a precision of ${\mathcal O}(\Delta t^{r'})$ then we can calculate $u_{k-1}$ with a precision of ${\mathcal O}(\Delta t^{r'+1})$ with the method $AB_{r'}$.

### Algorithm

- With the method $AB_1$, from $u_0$ we calculate $u_{-1}$ with a precision of ${\mathcal O}(\Delta t^2)$
- With the method $AB_2$, starting from $u_{0}$ and $u_{-1}$, we calculate $u_{1}$ with a precision of ${\mathcal O}(\Delta t^3)$
- For $r' = $3 to $r' = r$.
    - For $k=1$ to $k=r'-1$
        - With the method $AB_{r'-1}$, from $u_{1-k}, u_{2-k}, \ldots,u_{r'-1-k}$, we calculate $u_{-k}$ with a precision of ${\mathcal O}(\Delta t^{r'})$
    - For $k=1$ to $k=r'-1$
         - With the method $AB_{r'}$, from $u_{k-1}, u_{k-2}, \ldots,u_{k-r'}$, we calculate $u_{k}$ with a precision of ${\mathcal O}(\Delta t^{r'+1})$

At the end of this algorithm, the values $u_0, u_1, \ldots u_{r-1}$ are known with a precision of ${\mathcal O}(\Delta t^{r+1})$, we can launch the algorithm $AB_r$.


### The Adams-Bashforth Method
For an Adams-Bashforth method of order $r$ in time and spectral $\tau$, we first introduce a mesh in the $\tau$ direction.\
$\tau_{\ell} = \ell \Delta \tau, \ell = 0, \ldots, N_{\tau}-1$. Where $N_{\tau}$ is the number of points of discretization. If we apply the Fourier transform to the two-scale equation, we obtain 

$$\frac{\partial \hat{U}_\ell}{\partial t} + \frac{i\ell}{\varepsilon}\hat{U}_\ell = \hat{F}_\ell(t), \;\; \ell=-N_\tau/2, \dots, N_\tau/2-1,$$

with

$$U(t, \tau_k) = \sum_{\ell=-N_{\tau}/2}^{N_{\tau}/2-1} \hat{U}_{\ell}(t) e^{i\ell k 2\pi/N_{\tau}} \;\;\; \text{ and } F(\tau_k, U(t, \tau_k)) = \sum_{\ell=-N_{\tau}/2}^{N_{\tau}/2-1} \hat{F}_{\ell}(t) e^{i\ell k 2\pi/N_\tau}.$$

and the inverse (discrete) Fourier transform formulae

$$\hat{U}_{\ell}(t) = \frac{1}{N_{\tau}}\sum_{k=0}^{N_{\tau}-1} U(t, \tau_k) e^{-i\ell k 2\pi/N_{\tau}}\;\;\; \text{ and } \hat{F}_{\ell}(t) = \frac{1}{N_{\tau}}\sum_{k=0}^{N_{\tau}-1} F(\tau_k, U(t, \tau_k))e^{-i\ell k 2\pi/N_{\tau}}.$$

If we wish to calculate $\hat{F}_{\ell}$ from $\hat{U}_{\ell}$ we have the following formula 

$$F(\tau_k,U(t, \tau_k)) = e^{-\tau_k A}f(e^{\tau_k A}U(t, \tau_k)) \;\;$$

from which the Fourier transform is calculated in $\tau$ from the discrete Fourier transform formulas above.

Now, given a time step $\Delta t>0$ and a discretization in time $t_n=n\Delta t$, we can write the following Duhamel formula ($n\geq 0$)

$$\hat{U}_{\ell}(t_{n+1}) 
= e^{-i\ell\Delta t/\varepsilon}\hat{U}_{\ell}(t_{n})  + \int_0^{\Delta t} e^{-i\ell(\Delta t -s)/\varepsilon} \hat{F}_\ell(t_n+s)ds.$$

Thus, to obtain a $(r+1)$ scheme, we can 
approaches the function $\hat{F}_\ell(t_n+s)$ by the Lagrange polynomial of order $r$ interpolator at points $t_{n-j}, j=0, \dots, r$. This polynomial is written 

$$\hat{F}_\ell(t_n+s) \approx \sum_{k=0}^r \Big(\Pi_{j=0, j\neq k}^r \frac{s+j \Delta t}{(j-k)\Delta t} \Big) \hat{F}_\ell(t_n-t_j), \;\; n\geq 0.$$

Thus, from this approximation, we integrate exactly, which requires the following formulas

$$p^{[r]}_{\ell, j} = \int_0^{\Delta t}e^{-i\ell(\Delta t -s)/\varepsilon}\Big( \Pi_{j=0, j\neq k}^r \frac{s+j \Delta t}{(j-k)\Delta t}\Big) ds,$$

for each $j$ and $\ell$ such that $0\leq j\leq r, \; \ell=-N_\tau/2, \dots, N_\tau/2-1$. These coefficients $p^{[r]}_{\ell, j}$ can be pre-calculated and stored once and for all. Thus, the schema is finally written

$$\hat{U}_{\ell}^{n+1}= e^{-i\ell\Delta t/\varepsilon}\hat{U}_{\ell}^n + \sum_{j=0}^r p^{[r]}_{\ell, j} \hat{F}_\ell^{n-j},$$

with $\hat{U}_{\ell}^n \approx \hat{U}_{\ell}(t_n)$ and $\hat{F}_\ell^{n-j}\approx \hat{F}_\ell(t_{n-j})$. 

We can verify that the truncation error in this schema is  ${\mathcal O}(\Delta t^{r+1})$, once the initial values  $\hat{U}_\ell^1, \dots, \hat{U}_\ell^r$ have been calculated.  

## Non-homogeneous case $f(u, t)$

Here we consider the case where $f$ depends on the variable $t$.

$$\frac{d u(t)}{dt} = \frac{1}{\varepsilon} A u(t) + f(u(t), t), \;\;\; u(t=t_{0})=u_{0}, \;\; \varepsilon\in ]0, 1] \;\;\;\; (3)$$

The non-homogeneous case (3) falls under (1), by entering the variable 
$\theta : t \in [t_{0}, t_{\text{end}}] \mapsto \theta(t) = t\in \mathbb{R}$ 
which allows us to reformulate the non-homogeneous case (3) into a homogeneous problem of the form (1).
Indeed, we rephrase (3) as follows

$$\begin{aligned}
\frac{d u(t) }{dt}  & = \frac{1}{\varepsilon} Au(t) + f(u(t), \theta(t)), \\
\frac{d \theta(t) }{dt}  & = 1
\end{aligned}\;\;\;\;(4)$$

with the initial condition  $u(t_{0})=u_{0}, \theta(t_{0})=t_{0}$. Thus, the problem (4) is rewritten into an equation\
satisfied by $y: t\in [t_{0}, t_{fin}] \mapsto y(t) =(u(t), \theta(t))\in \mathbb{R}^{n+1}$

$$\frac{d y}{dt} = \frac{1}{\varepsilon} \tilde{A} y + g(y), \;\; y(t_{0})=(u_{0}, t_{0}),$$

with $\tilde{A}\in{\mathcal M}_{n+1, n+1}(\mathbb{R})$

$$\tilde{A}=
\left(
\begin{array}{cccc}
   &    &     & 0 \\
   & A &     & 0 \\
   &    &    & 0 \\
0 & 0 & 0 &  0 
\end{array}
\right) \;\;\;\; \text{ and } \;\;\;\;
g(y)=g(u, \theta) = \left(
\begin{array}{cccccc}
f(u, \theta) \\
1
\end{array}
\right) \in \mathbb{R}^{n+1}.$$

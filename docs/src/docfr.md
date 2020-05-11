

L'objectif de ce package Julia consiste à valoriser les développements récents effectués au sein de MINGuS sur les méthodes numériques Uniformément Précises (UA pour Uniformly Accurate) pour les problèmes hautement oscillants. On se propose de résoudre l'équation suivante 
$$
\frac{d u(t)}{dt} = \frac{1}{\varepsilon} A u(t) + f(u(t)), \;\;\; u(t=t_{0})=u_{0}, \;\; \varepsilon\in ]0, 1], \;\;\;\;\;\;\;\;\;\;\;\; (1)
$$
avec 
- $u : t\in [t_{0}, t_{fin}] \mapsto u(t)\in \mathbb{R}^n$, $t_{0}, t_{fin}\in \mathbb{R}$, 
- $u_{0}\in \mathbb{R}^n$, 
- $A\in {\cal M}_{n,n}(\mathbb{R})$ est telle que $\tau \mapsto \exp(\tau A)$ est périodique,  
- $f : u\in  \mathbb{R}^n \mapsto \mathbb{R}^n$.

Il s'agit ici d'écrire une documentation explicative du package *Julia* contenant la méthode double-échelle (voir [^1], [^2] et [^3]). Ce package s'inspire du package Differential Equations  qui est considéré comme un des meilleurs packages *Julia* existant. 


# Méthode numérique
## Formulation double-échelle 
Dans un premier temps, on réécrit l'équation (1) à l'aide du  changement d'inconnue $w(t)=\exp(-(t-t_{0})A/\varepsilon) u(t)$ pour obtenir 
$$
\frac{d w(t)}{dt} = F\Big( \frac{t-t_{0}}{\varepsilon}, w(t) \Big), \;\;\; w(t_{0})=u_{0}, \;\; \varepsilon\in ]0, 1],
$$
où la fonction $F$ s'exprime à partir des données du problème original (1) 
$$
F\Big( \frac{s}{\varepsilon}, w \Big) = \exp(-sA/\varepsilon) \; f( \exp(sA/\varepsilon) \; w). 
$$
On introduit alors la fonction $U(t, \tau), \tau\in [0, 2\pi]$ telle que $U(t, \tau=(t-t_{0})/\varepsilon) = w(t)$. La fonnction double-échelle est alors solution de l'équation suivante 
$$
\frac{\partial U}{\partial t} + \frac{1}{\varepsilon} \frac{\partial U}{\partial \tau} =  F( \tau, U), \;\;\; U(t=t_{0}, \tau)=\Phi(\tau), \;\; \varepsilon\in ]0, 1], \;\;\;\;\;\;\;\;\;\; (2)
%\label{eq_tau}
$$
où $\Phi$ est une fonction vérifiant $\Phi(\tau=0)=u_{0}$ choisie de sorte que la solution $U$ de (2) est régulière (voir [^1] et [^2]). 


## Discrétisation 
La méthode numérique est basée sur une discrétisation de l'équation (2). Dans la direction $\tau$, une méthode spectrale est utilisée alors que pour le temps $t$, une méthode de type exponentielle Adams-Bashforth permet de construire une méthode d'ordre élevé (voir [^2]). L'initialisation repose sur une technique "butterfly" (aller-retour autour du temps initial). 

### Initialisation 
Soit r l'ordre de la méthode $AB_r$.
Soit $\Delta t$ le pas de temps, pour $i \in \{ -r, -(r-1), \ldots, r-1, r\}$,
nous notons $u_i = u(t_0+i \Delta t)$.
Soit $r'$ les ordres des méthodes AB intermédiaires que nous allons utiliser.
Si $u_{k}$ est connue avec une précision de ${\cal O}(\Delta t^{r'+1})$, et pour $r' \geq 2, u_{k-1}, \ldots, u_{k-r'+1}$ sont connues avec une précision ${\cal O}(\Delta t^{r'})$ alors on peut calculer $u_{k+1}$ avec une précision de ${\cal O}(\Delta t^{r'+1})$ avec la méthode $AB_{r'}$
De même, si $u_{k}$ est connue avec une précision de ${\cal O}(\Delta t^{r'+1})$, et pour $r' \geq 2, u_{k+1}, \ldots, u_{k+r'-1}$ sont connues avec une précision ${\cal O}(\Delta t^{r'})$ alors on peut calculer $u_{k-1}$ avec une précision de ${\cal O}(\Delta t^{r'+1})$ avec la méthode $AB_{r'}$.
<!--- #### Algorithme(old)
Notons que, dans l'algorithme qui suit, l'ordre 2 est détaillé pour faciliter la compréhension, cet ordre aurait pu entrer dans le cadre général.
- Avec la méthode $AB_1$, à partir de $u_0$ on calcule $u_1$ avec une précision de ${\cal O}(\Delta t^2)$
- Avec la méthode $AB_2$, à partir de $u_0$ et $u_1$, on calcule $u_{-1}$ avec une précision de ${\cal O}(\Delta t^3)$
- Avec la méthode $AB_2$, à partir de $u_0$ et $u_{-1}$, on calcule $u_{1}$ avec une précision de ${\cal O}(\Delta t^3)$
- Avec la méthode $AB_2$,  à partir de $u_1$ et $u_0$, on calcule $u_{2}$ avec une précision de ${\cal O}(\Delta t^3)$
- Pour $r' = 3$ à $r' = r$
    - Pour $k=1$ à $k=r'-1$
        - Avec la méthode $AB_{r'}$, à partir de $u_{1-k}, u_{2-k}, \ldots,u_{r'-k}$, on calcule $u_{-k}$ avec une précision de ${\cal O}(\Delta t^{r'+1})$
    - Pour $k=1$ à $k=r'$
         - Avec la méthode $AB_{r'}$, à partir de $u_{k-1}, u_{k-2}, \ldots,u_{k-r}$, on calcule $u_{k}$ avec une précision de ${\cal O}(\Delta t^{r'+1})$

A la fin de cet algorithme, les valeurs $u_1, \ldots u_r$ sont connues avec une précision de ${\cal O}(\Delta t^{r+1})$, nous avons même une valeur de plus que nécessaire pour lancer l'algorithme $AB_r$.
-->
### Algorithme
- Avec la méthode $AB_1$, à partir de $u_0$ on calcule $u_{-1}$ avec une précision de ${\cal O}(\Delta t^2)$
- Avec la méthode $AB_2$, à partir de $u_{0}$ et $u_{-1}$, on calcule $u_{1}$ avec une précision de ${\cal O}(\Delta t^3)$
- Pour $r' = 3$ à $r' = r$
    - Pour $k=1$ à $k=r'-1$
        - Avec la méthode $AB_{r'-1}$, à partir de $u_{1-k}, u_{2-k}, \ldots,u_{r'-1-k}$, on calcule $u_{-k}$ avec une précision de ${\cal O}(\Delta t^{r'})$
    - Pour $k=1$ à $k=r'-1$
         - Avec la méthode $AB_{r'}$, à partir de $u_{k-1}, u_{k-2}, \ldots,u_{k-r'}$, on calcule $u_{k}$ avec une précision de ${\cal O}(\Delta t^{r'+1})$

A la fin de cet algorithme, les valeurs $u_0, u_1, \ldots u_{r-1}$ sont connues avec une précision de ${\cal O}(\Delta t^{r+1})$, nous pouvons lancer l'algorithme $AB_r$.
   


### La méthode Adams-Bashforth
Pour une méthode d'Adams-Bashforth d'ordre $r$ en temps et spectral en $\tau$, on introduit tout d'abord un maillage dans la direction $\tau$ 
$\tau_\ell = \ell \Delta \tau, \; \ell=0, \dots, N_\tau$ avec $\Delta\tau=2\pi/N_\tau$, $N_\tau$ étant le nombre de points de discrétisation. Si on applique la transformée de Fourier à l'équation double-échelle, on obtient 
$$
\frac{\partial \hat{U}_\ell}{\partial t} + \frac{i\ell}{\varepsilon}\hat{U}_\ell = \hat{F}_\ell(t), \;\; \ell=-N_\tau/2, \dots, N_\tau/2-1, 
$$
avec 
$$
U(t, \tau) = \sum_{\ell=-N_\tau/2}^{N_\tau/2-1} \hat{U}_\ell(t) e^{i\ell\tau} \;\;\; \mbox{ et } F(\tau, U(t, \tau)) = \sum_{\ell=-N_\tau/2}^{N_\tau/2-1} \hat{F}_\ell(t) e^{i\ell\tau}. 
$$
Si on souhaite calculer $\hat{F}$ à partir de $\hat{U}$ nous avons la formule <font color=red> à vérifier</font>
$$
\hat{F}_{\ell}(t) = \sum_{r=-N_\tau/2}^{N_\tau/2-1} e^{- \frac{\pi i r \ell}{N_\tau}}  e^{\frac{\ell}{\varepsilon} A} F( t, \frac{e^{-\frac{\ell}{\varepsilon} A}}{N_\tau}\sum_{k=-N_\tau/2}^{N_\tau/2-1} e^{\frac{\pi i k \ell}{N_\tau}} \hat{U}_k(t))
$$

Maintenant, étant donné un pas de temps $\Delta t>0$ et une discrétisation en temps $t_n=n\Delta t$, on peut écrire la formule de Duhamel suivante ($n\geq 0$) 
$$
\hat{U}_{\ell}(t_{n+1}) 
= e^{-i\ell\Delta t/\varepsilon}\hat{U}_{\ell}(t_{n})  + \int_0^{\Delta t} e^{-i\ell(\Delta t -s)/\varepsilon} \hat{F}_\ell(t_n+s)ds. 
$$
Ainsi, pour obtenir un schéma d'ordre $(r+1)$, on 
approche la fonction $\hat{F}_\ell(t_n+s)$ par le polynôme de Lagrange d'ordre $r$ interpolateur aux points $t_{n-j}, j=0, \dots, r$. Ce polynôme s'écrit 
$$
\hat{F}_\ell(t_n+s) \approx \sum_{k=0}^r \Big(\Pi_{j=0, j\neq k}^r \frac{s+j \Delta t}{(j-k)\Delta t} \Big) \hat{F}_\ell(t_n-t_j), \;\; n\geq 0.  
$$
Ainsi, à partir de cette approximation, on intègre exactement, ce qui demande les formules suivantes 
$$
p^{[r]}_{\ell, j} = \int_0^{\Delta t}e^{-i\ell(\Delta t -s)/\varepsilon}\Big( \Pi_{j=0, j\neq k}^r \frac{s+j \Delta t}{(j-k)\Delta t}\Big) ds,   
$$
pour chaque $j$ et $\ell$ tels que $0\leq j\leq r, \; \ell=-N_\tau/2, \dots, N_\tau/2-1$. Ces coefficients $p^{[r]}_{\ell, j}$ peuvent être pré-calculés et stockés une fois pour toute. Ainsi, le schéma s'écrit finalement 
$$
\hat{U}_{\ell}^{n+1}= e^{-i\ell\Delta t/\varepsilon}\hat{U}_{\ell}^n + \sum_{j=0}^r p^{[r]}_{\ell, j} \hat{F}_\ell^{n-j},  
$$
avec $\hat{U}_{\ell}^n \approx \hat{U}_{\ell}(t_n)$ et $\hat{F}_\ell^{n-j}\approx \hat{F}_\ell(t_{n-j})$. 

On peut vérifier que l'erreur de troncature de ce schéma est ${\cal O}(\Delta t^{r+1})$, une fois que les valeurs initiales $\hat{U}_\ell^1, \dots, \hat{U}_\ell^r$ ont été calculées.  



## Cas non-homogène $f(u, t)$
On considère ici le cas où $f$ dépend de la variable $t$

$$
\frac{d u(t)}{dt} = \frac{1}{\varepsilon} A u(t) + f(u(t), t), \;\;\; u(t=t_{0})=u_{0}, \;\; \varepsilon\in ]0, 1], 
%\tag{eq1_t}
$$

Le cas non-homogène \ref{eq1_t} entre dans le cadre \eqref{eq1} 
en introduisant la variable $\theta : t \in [t_{0}, t_{fin}] \mapsto \theta(t) = t\in \mathbb{R}$ 
qui permet de reformuler le cas non-homogène \eqref{eq1_t} en un problème homogène de la forme \eqref{eq1}. 
En effet, on reformule \eqref{eq1_t} de la manière suivante 
$$
\begin{aligned}
\frac{d u(t) }{dt}  & = \frac{1}{\varepsilon} Au(t) + f(u(t), \theta(t)), \\
\frac{d \theta(t) }{dt}  & = 1,
\end{aligned}
$$
%\label{eq1_t_theta}

avec la condition initiale $u(t_{0})=u_{0}, \theta(t_{0})=t_{0}$. Ainsi, le problème \eqref{eq1_t_theta} se réécrit en une équation 
satisfaite par $y: t\in [t_{0}, t_{fin}] \mapsto y(t) =(u(t), \theta(t))\in \mathbb{R}^{n+1}$
$$
\frac{d y}{dt} = \frac{1}{\varepsilon} \tilde{A} y + g(y), \;\; y(t_{0})=(u_{0}, t_{0}), 
$$
avec $\tilde{A}\in{\cal M}_{n+1, n+1}(\mathbb{R})$
$$
\tilde{A}=
\left(
\begin{array}{cccc}
   &    &     & 0 \\
   & A &     & 0 \\
   &    &    & 0 \\
0 & 0 & 0 &  0 
\end{array}
\right) \;\;\;\; \mbox{ et } \;\;\;\;
g(y)=g(u, \theta) = \left(
\begin{array}{cccccc}
f(u, \theta) \\
1
\end{array}
\right) \in \mathbb{R}^{n+1}. 
$$

## La méthode Runge-Kutta "classique" (ordre 4) adaptée à l'exponentielle

Notations : 
- Nous notons $G$ la fonction qui permet de passer de $\hat{U}$ à $\hat{f}$, ainsi $\hat{f} = G(\hat{U})$.
- Nous notons $S_{t_0}^{t_1}(t_2,\ell)$ l'intégrale $S_{t_0}^{t_1}(t_2,\ell) = \int_{t_0}^{t_1} e^{- i \ell (t_2 - s)/\varepsilon} ds = ( i \varepsilon / \ell) ( e^{- i \ell (t_2 - t_1)/\varepsilon}-e^{- i \ell (t_2 - t_0)/\varepsilon})$


Voici les calculs


- $u_{1,\ell} = \hat{U}_{n, \ell}$
- $u_{2,\ell} = e^{- i \ell h_n /(2 \varepsilon)}\hat{U}_{n, \ell} + S_0^{h_n /2} ( h_n /2,\ell ) G_{\ell}(u_1)$
- $u_{3,\ell} = e^{- i \ell h_n /(2 \varepsilon)}\hat{U}_{n, \ell} + S_0^{h_n /2} ( h_n /2,\ell )  G_{\ell}(u_2)$
- D'après (28) du papier, avec $c=-i \ell h_n /\varepsilon$, on a   
- $u_{4,\ell} = e^{- i \ell h_n /(2\varepsilon)}u_{2,\ell} + S_0^{h_n/2} ( h_n/2,\ell )[ 2 G_{\ell}(u_3)-G_{\ell}(u_1)]$

$$\hat{U}_{n+1, \ell} = e^{- i \ell h_n /\varepsilon}\hat{U}_{n, \ell} + S_0^{h_n/6} ( h_n,\ell ) G_{\ell}(u_1) + \frac{1}{2} S_{h_n/6}^{5 h_n/6} ( h_n,\ell ) ( G_{\ell}(u_2) + G_{\ell}(u_3) )\\ + S_{5h_n/6}^{ h_n} ( h_n,\ell ) G_{\ell}(u_4)$$

D'après (29) du papier, avec $c=-i \ell h_n /\varepsilon$, on a
$$
\hat{U}_{n+1, \ell} = e^{- i \ell h_n /\varepsilon}\hat{U}_{n, \ell} +  G_{\ell}(u_1) [-4+i \ell h_n /\varepsilon + e^{-i \ell h_n /\varepsilon}(4+3i \ell h_n /\varepsilon+(i \ell h_n /\varepsilon)^2]\\+ (2 G_{\ell}(u_2) + G_{\ell}(u_3) )[-2-i \ell h_n /\varepsilon+e^{-i \ell h_n /\varepsilon}(2-i \ell h_n /\varepsilon)]\\ + G_{\ell}(u_4)[-4+3i \ell h_n /\varepsilon -(i \ell h_n /\varepsilon)^2 + e^{-i \ell h_n /\varepsilon}(4+i \ell h_n /\varepsilon)]/(h_n^2 (i \ell h_n /\varepsilon)^3)
$$


# Utilisation



## Paramètres d'entrée
Les arguments d'entrée utilisent le même format que le package ODE. 

Ainsi, dans un premier temps, il faut définir les arguments nécessaires pour construire le problème  \eqref{eq1}, à savoir

- la fonction $f$ (sous la forme *Julia*) 
- la condition initiale $u_{0}$
- les temps initial $t_{0}$ et final $t_{fin}$ 
- le deuxième paramètre de la fonction 
- la matrice $A$ 
- $\varepsilon \in ]0, 1]$ 

Note : *Vous devez saisir "à la main" le ] et ^C, pour le reste vous pouvez copier-coller*.

```jl
]
add https://gitlab.inria.fr/ua/HighlyOscillatoryProblems.jl.git
^C
using HighlyOscillatoryProblems
A=[0 0 1 0 ; 0 0 0 0 ; -1 0 0 0 ; 0 0 0 0]
fct = (u,p,t) ->  [ 0, u[4], 2*u[1]*u[2], -u[2] - u[1]^2 + u[2]^2 ] 
epsilon= 0.0001
t_min=0.0
t_max=3.0
u0 = [0.55, 0.12, 0.03, 0.89]
prob = HiOscODEProblem(fct, u0, (t_min,t_max), missing, A, epsilon) 
```

A partir du problème ```prob```, on peut maintenant passer à sa résolution numérique. 
Pour cela, on définit les paramètres numériques 
- le nombre de tranches de temps $N_t$ qui définit le pas de temps $\Delta t = \frac{t_0-t_{fin}}{N_t}$ 
- l'ordre $r$ de la méthode 
- le nombre de points $N_\tau$ dans la direction $\tau$ 
- l'ordre de la préparation $q$ de la condition initiale 
- précision  <font color=red> non encore implémentée</font>. 

Par défaut, les paramètres sont : $N_t=100$, $r=4$, $N_\tau=32$ et $q=r+2=6$
Pour résoudre le problème avec les paramètres par défaut, il suffit d'appeler la commande `solve` avec le problème déjà défini comme paramètre
```jl     
sol = solve(prob) 
```
Ce qui est équivalent à cet appel 
```jl     
sol = solve(
    prob;
    nb_tau=32, 
    order=4, 
    order_prep=order+2,
    dense=true, 
    nb_t=100, 
    getprecision=dense, 
    verbose=100,    
    par_u0=missing,
    p_coef=missing,
) 
```
### Définition exhaustive des paramètres
- `prob` : problème défini par `HiOscODEProblem` 
- `nb_tau=32` : $N_{\tau}$
- `order=4` : ordre $r$ de la méthode
- `order_prep=order+2` : ordre de préparation des données initiales
- `dense=true` : indique si l'on conserve ou non les données issues de la transformée de fourier, si `dense=false`, le traitement est plus rapide mais l'interpolation ne peut plus être faite.
- `nb_t=100` : $N_t$
- `getprecision=dense` : indique si la précision est calculée, la méthode utilisée pour calculer la précision multiplie par 2 le temps de traitement.
- `verbose=100` : niveau de trace, si `verbose=0` alors rien n'est affiché.
- `par_u0` : données initiales préparées, si on doit faire plusieurs appels à `solve` avec les mêmes données initiales et au même ordre on peut passer en paramêtre les données déjà calculées.
- `p_coef` : tableau avec les coefficients de la méthode Adams-Bashforth. Ce tableau peut être utilisé pour optimiser plusieurs appels avec les mêmes paramètres. 
## Arguments de sortie
En sortie, une structure de type `HiOscODESolution`.
Cette structure peut être vue comme une fonction de t, elle peut aussi être vue comme un tableau de taille $N_t + 1$. Cette structure contient aussi les champs `absprec` et `relprec` qui sont les précisions, respectivement absolue et relative, calculées.
### Exemple
```jl     
julia> sol = solve(prob);
solve function prob=HiOscODEProblem with uType Array{Float64,1} and tType Float64. In-place: nothing
timespan: (0.0, 3.0)
u0: [0.55, 0.12, 0.03, 0.89],
 nb_tau=32, order=4, order_prep=6, dense=true,
 nb_t=100, getprecision=true, verbose=100

x 100/100

 99/99

julia> t=2.541451547
2.541451547

julia> sol(t)
4-element Array{Float64,1}:
 -0.536667845897295
  1.593257176840297
 -0.12420061944907212
  0.7184374612958457

julia> sol[end]
4-element Array{Float64,1}:
  0.36316109321808354
  2.0379196624858955
 -0.4141248638226731
  1.3087136174628513

julia> sol(3.0)
4-element Array{Float64,1}:
  0.36316109321808626
  2.037919662485913
 -0.4141248638226772
  1.308713617462862

julia> sol.absprec
2.4721022528746903e-5

julia> sol.relprec
9.952927361881597e-6


```
Pour visualiser le résultat, on peut aussi utiliser Plot, par exemple
```jl     
using Plots
plot(sol) 
```
Ce qui donne
![](https://codimd.math.cnrs.fr/uploads/upload_a9389e1a0443b2e838dbfc7fb8658cc0.png)
Si l'on souhaite sauvegarder le résultat dans un fichier (pdf ou png) voici les commandes
```jl     
using Plots
p = plot(sol)
savefig(p,"out/plot.png")
```

# Exemples 
## Hénon-Heiles
On considère le système de Hénon-Heiles satisfait par $u(t)=(u_1, u_2, u_3, u_4)(t)$ 
$$
\frac{d u }{dt} = \frac{1}{\varepsilon} Au + f(u), \;\;\; u(0)=u_0\in\mathbb{R}^4, 
$$
où $A$ et $f$ sont choisis comme suit 
$$
A=
\left(
\begin{array}{cccc}
0 & 0 & 1 & 0  \\
0 & 0 & 0 & 0  \\
-1 & 0 & 0 & 0  \\
0 & 0 & 0 & 0  
\end{array}
\right) \;\;\;\; \mbox{ et } \;\;\;\;
f(u) = \left(
\begin{array}{cccc}
0 \\
u_4\\
-2 u_1 u_2\\
-u_2-u_1^2+u_2^2
\end{array}
\right). 
$$
on choisit par exemple, $\varepsilon=0.001$
et $u_0 = (0.12, 0.12, 0.12, 0.12)$
```jl
using HighlyOscillatoryProblems
A=[0 0 1 0 ; 0 0 0 0 ; -1 0 0 0 ; 0 0 0 0]
fct = (u,p,t) ->  [ 0, u[4], 2*u[1]*u[2], -u[2] - u[1]^2 + u[2]^2 ] 
epsilon= 0.001
t_min=0.0
t_max=1.0
u0 = [0.12, 0.12, 0.12, 0.12]
prob = HiOscODEProblem(fct, u0, (t_min,t_max), missing, A, epsilon)
sol = solve(prob);
using Plots
plot(sol)
```

![](https://codimd.math.cnrs.fr/uploads/upload_1514c7513ad489a27a3a3287dfcb2666.png)


## Particule chargée

On considère un système de particules chargées sous l'effet d'un champ électro-magnétique extérieur 
$(E(t, x), B(t, x))\in \mathbb{R}^6$. Les particules sont décrites de manière dynamique par leur position 
$x(t)\in\mathbb{R}^3$ et leur vitesse $v(t)\in\mathbb{R}^3$. On indexera par  $i$  la $i$-ème composante d'un vecteur.
Les équations de Newton appliquées à une particule s'écrivent  
\begin{eqnarray*}
\frac{d x(t) }{dt}&=& v(t) \\
\frac{d v(t) }{dt}&=& \frac{e}{m} \left[E(t, x(t)) + v(t)\times B(t, x(t))\right]. 
\end{eqnarray*}
On va supposer que le champ magnétique s'écrit $B(t, x)=(0, 0, 1)^T$ et sous un certain scaling, on considère l'équation 
suivante 
\begin{eqnarray*}
\frac{d x_1(t) }{dt} &=& \frac{1}{\varepsilon}v_1(t) \\
\frac{d x_2(t) }{dt} &=& \frac{1}{\varepsilon} v_2(t) \\
\frac{d x_3(t) }{dt}&=& v_3(t) \\
\frac{d v_1(t) }{dt} &=& E_1(t, x(t)) + \frac{1}{\varepsilon}v_2(t)\\
\frac{d v_2(t) }{dt} &=& E_2(t, x(t)) - \frac{1}{\varepsilon}v_1(t)\\
\frac{d v_3(t) }{dt} &=& E_3(t, x(t)) 
\end{eqnarray*}
ce qui se réécrit sous la forme suivante 
$$
\frac{d u(t) }{dt}= \frac{1}{\varepsilon}A u(t) + F(t, u(t)), 
$$
où le vecteur inconnu $u(t)=(x(t), v(t))\in\mathbb{R}^6$, $A$ est une matrice carrée de taille $6\times 6$ 
et $F$ est une fonction à valeur dans $\mathbb{R}^6$. $A$ et $F$ sont données par 
$$
A=
\left(
\begin{array}{cccccc}
0 & 0 & 0 & 1 & 0 & 0 \\
0 & 0 & 0 & 0 & 1 & 0 \\
0 & 0 & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 & 1 & 0 \\
0 & 0 & 0 & -1 & 0 & 0 \\
0 & 0 & 0 & 0 & 0 & 0 
\end{array}
\right) \;\;\;\; \mbox{ et } \;\;\;\;
F(t, u(t)) = \left(
\begin{array}{cccccc}
0 \\
0\\
u_6(t)\\
E_1(t, u_1(t), u_2(t), u_3(t))\\
E_2(t,  u_1(t), u_2(t), u_3(t)\\
E_3(t,  u_1(t), u_2(t), u_3(t)
\end{array}
\right). 
$$
On peut considérer la fonction $E=(E_1, E_2, E_3)$ suivante 
$$
E(t, x) =
\left(
\begin{array}{ccc}
\cos(x_1/2)\sin(x_2)\sin(x_3)/2\\
\sin(x_1/2)\cos(x_2)\sin(x_3)\\
\sin(x_1/2)\sin(x_2)\cos(x_3)
\end{array}
\right)
$$
```jl
using HighlyOscillatoryProblems
A = [0 0 0 1 0 0; 0 0 0 0 1 0;0 0 0 0 0 0; 0 0 0 0 1 0; 0 0 0 -1 0 0; 0 0 0 0 0 0]
function fparticle(u, p, t)
    s1, c1 = sincos(u[1]/2)
    s2, c2 = sincos(u[2])
    s3, c3 = sincos(u[3])
    return [0, 0, u[6], c1*s2*s3/2, s1*c2*s3, s1*s2*c3]
end
epsilon= 0.05
t_min=0.0
t_max=1.0
u0 = [1.0, 1.5, -0.5, 0, -1.2, 0.8]
prob = HiOscODEProblem(fparticle, u0, (t_min,t_max), missing, A, epsilon)
sol=solve(prob)
plot(sol)
```
![](https://codimd.math.cnrs.fr/uploads/upload_9880503c3c87a8d26e31405fed411f17.png)

## Cas non-homogène linéaire 
On considère le système linéaire non-homogène suivant satisfait par $u(t)=(u_1, u_2, u_3, u_4)(t)$ 
$$
\frac{d u }{dt} = \frac{1}{\varepsilon} Au + f(t, u), \;\;\; u(0)=u_0\in\mathbb{R}^4, 
$$
où $A$ et $f$ sont choisis comme suit 
$$
A=
\left(
\begin{array}{cccc}
0 & 0 & 1 & 0  \\
0 & 0 & 0 & 0  \\
-1 & 0 & 0 & 0  \\
0 & 0 & 0 & 0  
\end{array}
\right) \;\;\; \mbox{ et } \;\;\;
f(t, u) = Bu +\alpha t +\beta \;\; \mbox{ avec  } \;\;
B\in {\cal M}_{4, 4}(\mathbb{R}), \alpha, \beta \in \mathbb{R}^4, 
$$
$B, \alpha, \beta$ sont choisis aléatoirement.
Nous souhaitons obtenir une grande précision, nous allons donc utiliser des réel de type BigFloat, il sont codés sur 256 bits par défaut ce qui donne une borne de précision d'environ $2^{-256} \approx 10^{-77}$. 
Nous comparons à la fin un résultat calculé avec un résultat exact. 
```jl
using HighlyOscillatoryProblems
using Random
Random.seed!(1111)
A=[0 0 1 0 ; 0 0 0 0 ; -1 0 0 0 ; 0 0 0 0]
B = 2rand(BigFloat, 4, 4) - ones(BigFloat, 4, 4)
alpha = 2rand(BigFloat, 4) - ones(BigFloat, 4)
beta = 2rand(BigFloat, 4) - ones(BigFloat, 4)
fct = (u,p,t)-> B*u + t*p[1] +p[2]
u0 = [big"0.5", big"-0.123", big"0.8", big"0.7"]
t_min=big"0.0"
t_max=big"1.0"
epsilon=big"0.017"
prob = HiOscODEProblem(fct, u0, (t_min,t_max), (alpha, beta), A, epsilon, B)
sol = solve(prob, nb_t=10000, order=8)
sol.absprec
t=big"0.9756534187771"
sol(t)-getexactsol(sol.par_u0.parphi, u0, t)
using Plots
Plots.plot(sol.t,sol.u_tr)
```
### Calcul de la solution exacte
Il s'agit de calculer la solution exacte $u(t)$ de l'équation suivante à l'instant $t$
$$
\frac{d u }{dt} = \frac{1}{\varepsilon} Au + Bu +\alpha t +\beta, \;\;\; u(0)=u_0\in\mathbb{R}^4\mbox{, } A \mbox{ et }B \mbox{ sont définies plus haut } 
$$
Soit $M = \frac{1}{\varepsilon} A + B$

$C = e^{-t_0 M}u_0 +M^{-1} e^{-t_0 M} (t_0\alpha+\beta)+ M^{-2} e^{-t_0 M} \alpha$
$C_t = -M^{-1} e^{-t M} (t\alpha+\beta)-M^{-2} e^{-t M} \alpha$

$u(t) = e^{t M} ( C + C_t)$

Ce qui, traduit en langage Julia, donne le code de la fonction `getexactsol` : 
```julia
function getexactsol(par::PreparePhi, u0, t)
    @assert !ismissing(par.matrix_B) "The debug matrix is not defined"
    sparse_A = par.sparse_Ap[1:(end-1),1:(end-1)]
    m = (1/par.epsilon)*sparse_A+par.matrix_B
    t0 = par.t_0
    if ismissing(par.paramfct)
        return exp((t-t0)*m)*u0
    end
    a, b = par.paramfct
    mm1 = m^(-1)
    mm2 = mm1^2
    e_t0 = exp(-t0*m)
    C = e_t0*u0 + mm1*e_t0*(t0*a+b)+mm2*e_t0*a
    e_inv = exp(-t*m)
    e = exp(t*m)
    C_t = -mm1*e_inv*(t*a+b)-mm2*e_inv*a
    return e*C+e*C_t
end

```


### Bibliographie
[^1]: P. Chartier, N. Crouseilles, M. Lemou, F. Méhats, 
Numer. Math., 129, pp. 211-250, (2015).

[^2]: P. Chartier, M. Lemou, F. Méhats, X. Zhao, 
submitted. 

[^3]: N. Crouseilles, M. Lemou, F. Méhats, 
J. Comput. Phys, 248, pp.287-308, (2013). 





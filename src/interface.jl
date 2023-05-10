using SciMLBase

using Reexport
@reexport using SciMLBase

include("twoscales_pure_ab.jl")

struct HOODEFunction{iip,N} <: SciMLBase.AbstractODEFunction{iip}
    f::Any
end
(fct::HOODEFunction{false,3})(u, p, t) = fct.f(u, p, t)
(fct::HOODEFunction{false,2})(u, p, t) = fct.f(u, p)
(fct::HOODEFunction{false,1})(u, p, t) = fct.f(u)
function (fct::HOODEFunction{true,4})(u, p, t)
    du = zero(u)
    fct.f(du, u, p, t)
    return du
end


"""
    HOODEProblem(f, u0, tspan, p, A, epsilon, B)

The HOODE problem is :

```math
\\frac{du}{dt} = ( \\frac{1}{\\epsilon} A + B)  u + f(u,p,t)
```

- The initial condition is ``u(tspan[1]) = u0``.
- The solution ``u(t)`` will be computed for ``tspan[1] ≤ t ≤ tspan[2]``
- Constant parameters to be supplied as the second argument of ``f``.
- Periodic Matrix of the problem.
- epsilon of the problem.
- Matrix of linear problem to get the exact solution

"""
struct HOODEProblem{T} <: SciMLBase.DEProblem
    f::HOODEFunction
    u0::Vector{T}
    tspan::Tuple{T,T}
    p::Any
    A::AbstractMatrix
    epsilon::T
    B::Union{Matrix,Missing}
    function HOODEProblem(
        f,
        u0::Vector{T},
        tspan::Tuple{T,T},
        p,
        A::AbstractMatrix,
        epsilon::T,
        B::Union{Matrix,Missing},
    ) where {T}
        fct = if hasmethod(f, (Vector{T}, Vector{T}, Any, T))
            HOODEFunction{true,4}(f)
        elseif hasmethod(f, (Vector{T}, Any, T))
            HOODEFunction{false,3}(f)
        elseif hasmethod(f, (Vector{T}, Any))
            HOODEFunction{false,2}(f)
        elseif hasmethod(f, (Vector{T},))
            HOODEFunction{false,1}(f)
        else
            println("err !!!!!")
        end
        return new{T}(fct, u0, tspan, p, A, epsilon, B)
    end # end of function
end # end of struct
function HOODEProblem(f, u0, tspan, p, A, epsilon)
    return HOODEProblem(f, u0, tspan, p, A, epsilon, missing)
end
struct HOODEInterpolation{T} <: SciMLBase.AbstractDiffEqInterpolation
    t::Vector{T}
    u_caret::Vector{Array{Complex{T},2}}
    parphi::PreparePhi
    order::Any
end
(interp::HOODEInterpolation)(t, idxs, deriv, p, continuity) = interp(t)
function (interp::HOODEInterpolation)(t)
    return _getresult(
        interp.t,
        interp.u_caret,
        t,
        interp.parphi,
        interp.t[1],
        interp.t[end],
        interp.order,
    )
    #    return _getresult(interp.u_caret, t, 
    #    interp.parphi, 
    #    interp.t[1], interp.t[end], 
    #    interp.order)
end
(interp::HOODEInterpolation)(vt::Vector{<:Number}) = RecursiveArrayTools.DiffEqArray(interp.(vt),vt)
if typeof(SciMLBase.AbstractTimeseriesSolution{Float64,Float64,Float64}) == DataType
    abstract type AbstractHOODESolution{T,N} <: SciMLBase.AbstractTimeseriesSolution{T,N,N} end
else
    abstract type AbstractHOODESolution{T,N} <: SciMLBase.AbstractTimeseriesSolution{T,N} end
end
struct HOODESolution{T} <: AbstractHOODESolution{T,T}
    u::Vector{Vector{T}}
    u_tr::Union{Vector{Vector{T}},Nothing}
    t::Vector{T}
    dense::Bool
    order::Integer
    par_u0::PrepareU0
    p_coef::CoefExpAB
    prob::HOODEProblem{T}
    retcode::Any
    interp::Union{HOODEInterpolation,Nothing}
    tslocation::Any
    absprec::Any
    relprec::Any
end
function HOODESolution(retcode::Symbol)
    return HOODESolution(
        nothing,
        nothing,
        nothing,
        nothing,
        nothing,
        nothing,
        retcode,
        nothing,
        nothing,
        nothing,
    )
end
(sol::HOODESolution)(t) = sol.dense ? sol.interp(t) : undef

abstract type AbstractHOODEAlgorithm end

struct HOODETwoScalesAB <: AbstractHOODEAlgorithm end
struct HOODEETDRK2 <: AbstractHOODEAlgorithm end
struct HOODEETDRK3 <: AbstractHOODEAlgorithm end
struct HOODEETDRK4 <: AbstractHOODEAlgorithm end

import SciMLBase:solve

"""
    function solve(prob::HOODEProblem{T}; 
    nb_tau::Integer=32, 
    order::Integer=4, 
    order_prep::Integer=order+2, 
    dense::Bool=true, 
    nb_t::Integer=100, 
    getprecision::Bool=dense,
    verbose=100,
    par_u0::Union{PrepareU0,Missing}=missing,
    p_coef::Union{CoefExpAB,Missing}=missing
    ) where T<:AbstractFloat

specific interface solver for Highly oscillatory problems, that an ODE of this form:
```math
\\frac{\\delta u(t)}{\\delta t} = \\frac{1}{\\varepsilon} A + F(u(t), t)
```
where ``u \\in \\R^n`` and  ``0 < \\varepsilon < 1``
``A`` must be a **periodic matrix** i.e. ``e^{t A} = e^{(t+\\pi) A}`` for any ``t \\in \\R``

## Argument :
- `prob::HOODEProblem{T}` : The problem to solve

## Keywords :
- `nb_tau::Integer=32` : number of values of FFT transform, must be power of twoscales_pure_ab
- `order::Integer=4` : order of Adams-Bashforth method, and also of the interpolatation
- `order_prep::Integer=order+2` : order of the preparation
- `dense::Bool=true` : if true it is possible to compute solution at any time of the interval
- `nb_t::Integer=100` : number of period slices
- `getprecision::Bool=dense` : compute the absolute and relative precision
- `par_u0::Union{PrepareU0,Missing}=missing` : preparation data for u0
- `p_coef::Union{CoefExpAB,Missing}=missing` : coefficients for Adams-Bashforth method

## Examples :
"""
function SciMLBase.solve(prob::HOODEProblem{T}; kwargs...) where {T<:AbstractFloat}
    return SciMLBase.solve(prob, HOODETwoScalesAB(); kwargs...)
end


function SciMLBase.solve(
    prob::HOODEProblem{T},
    alg::HOODETwoScalesAB;
    nb_tau::Integer = 32,
    order::Integer = 4,
    order_prep::Integer = order + 2,
    dense::Bool = true,
    nb_t::Integer = 100,
    getprecision::Bool = dense,
    verbose = 100,
    par_u0::Union{PrepareU0,Missing} = missing,
    p_coef::Union{CoefExpAB,Missing} = missing,
) where {T<:AbstractFloat}
    retcode = SciMLBase.ReturnCode.Success
    nb_tau = prevpow(2, nb_tau)
    traceln(
        100,
        "solve function prob=$prob,\n nb_tau=$nb_tau, order=$order, order_prep=$order_prep, dense=$dense,\n nb_t=$nb_t, getprecision=$getprecision, verbose=$verbose";
        verbose = verbose,
    )
    if getprecision
        s1 = SciMLBase.solve(
            prob;
            nb_tau = nb_tau,
            order = order,
            order_prep = order_prep,
            dense = dense,
            nb_t = nb_t,
            getprecision = false,
            verbose = verbose - 10,
            par_u0 = par_u0,
        )
        n_nb_t = Integer(floor(1.1373 * nb_t + 1))
        s2 = SciMLBase.solve(
            prob;
            nb_tau = nb_tau,
            order = order,
            order_prep = order_prep,
            dense = dense,
            nb_t = n_nb_t,
            getprecision = false,
            verbose = verbose - 10,
            par_u0 = s1.par_u0,
        )
        absprec = norm(s1.u[end] - s2.u[end])
        relprec = absprec / max(norm(s1.u[end]), norm(s2.u[end]))
        # if dense
        #     for i = 1:10
        #         t = rand(T) * (prob.tspan[2]-prob.tspan[1]) 
        #             + prob.tspan[1]
        #         a, b = s1(t), s2(t)
        #         ap = norm(a-b)
        #         rp = ap/max(norm(a),norm(b))
        #         absprec = max(absprec,ap)
        #         relprec = max(relprec,rp)
        #     end
        # end
        return HOODESolution(
            s1.u,
            s1.u_tr,
            s1.t,
            dense,
            order,
            s1.par_u0,
            s1.p_coef,
            prob,
            retcode,
            s1.interp,
            0,
            Float64(absprec),
            Float64(relprec),
        )
    end

    par_u0 = if ismissing(par_u0)
        parphi = PreparePhi(
            prob.epsilon,
            nb_tau,
            prob.A,
            prob.f,
            prob.B;
            t_0 = prob.tspan[1],
            paramfct = prob.p,
        )
        PrepareU0(parphi, order_prep, prob.u0)
    else
        par_u0
    end
    pargen = PrepareTwoScalesPureAB(
        nb_t,
        prob.tspan[2],
        order,
        par_u0,
        p_coef = p_coef,
        verbose = verbose,
    )
    if dense
        u_mat, _, u_caret = twoscales_pure_ab(pargen, res_fft = true, verbose = verbose)
        t = collect(0:nb_t) * (prob.tspan[2] - prob.tspan[1]) / nb_t .+ prob.tspan[1]
        interp = HOODEInterpolation{T}(t, u_caret, par_u0.parphi, order)
        u = reshape(mapslices(x -> [x], u_mat, dims = 1), size(u_mat, 2))
        u_tr = reshape(mapslices(x -> [x], transpose(u_mat), dims = 1), size(u_mat, 1))
    else
        u = Array{Array{T,1},1}(undef, 2)
        u[1] = copy(prob.u0)
        u[end] = twoscales_pure_ab(pargen, verbose = verbose, only_end = true)
        t = [prob.tspan[1], prob.tspan[2]]
        u_tr = nothing
        interp = nothing
    end
    return HOODESolution(
        u,
        u_tr,
        t,
        dense,
        order,
        par_u0,
        pargen.p_coef,
        prob,
        retcode,
        interp,
        0,
        nothing,
        nothing,
    )
end

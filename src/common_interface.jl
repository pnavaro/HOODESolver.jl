using LinearAlgebra
using GenericSchur

include("interface.jl")

struct LinearHOODEOperator{T} <: SciMLBase.AbstractDiffEqLinearOperator{T}
    epsilon::T
    A::AbstractMatrix
end
SciMLBase.isinplace(linop::LinearHOODEOperator, n::Int) = false

function LinearHOODEOperator(mat::AbstractMatrix{T}) where {T}
    P = eigvecs(mat .+ 0im) # ca bug avec les reels
    Ad = inv(P) * mat * P
    epsilon = 1 / maximum(abs.(Ad))
    A = epsilon * mat
    Aint = round.(A)
    @assert norm(A - Aint, Inf) <= eps(T) "The A matrix must be with integers A=$A"
    epsilon = norm(Aint, Inf) / norm(mat, Inf)
    LinearHOODEOperator{T}(epsilon, Aint)
end

*(v::T, L::LinearHOODEOperator{T}) where {T<:Number} =
    LinearHOODEOperator(L.epsilon / v, L.A)
Base.@propagate_inbounds Base.convert(::Type{AbstractMatrix}, L::LinearHOODEOperator) =
    (1 / L.epsilon) * L.A

import SciMLBase: isconstant

isconstant(_::LinearHOODEOperator) = true

function LinearHOODEOperator(linop::DiffEqArrayOperator)
    linop.update_func != SciMLBase.DEFAULT_UPDATE_FUNC &&
        error("no update operator function for HOODEAB Alg")
    LinearHOODEOperator(linop.A)
end

LinearHOODEOperator(linop::LinearHOODEOperator) = linop

function LinearHOODEOperator(odefct::ODEFunction)
    return LinearHOODEOperator(odefct.f)
end

isSplitODEProblem(probtype::Any) = false
isSplitODEProblem(probtype::SplitODEProblem) = true


"""
    HOODEAB( order=4, ntau=32)

Algorithm for High-Oscillatory equation
"""
struct HOODEAB{order,ntau} <: SciMLBase.AbstractODEAlgorithm
    HOODEAB(order::Int = 4; ntau = 32) = new{order,ntau}()
end
export HOODEAB

# OtherHOODE : structure used to store HOODESolution data in the dstats fields
struct OtherHOODE
    par_u0::PrepareU0
    p_coef::CoefExpAB
    absprec::Union{Nothing,Float64}
    relprec::Union{Nothing,Float64}
end

"""
    solve(prob::ODEProblem, alg::HOODEAB{order, ntau}; 
    dt=nothing,
    kwargs...
    ) where {order,ntau}

common interface solver for Highly oscillatory problems, that an ODE of this form

```math
\\frac{\\delta u(t)}{\\delta t} = \\frac{1}{\\varepsilon} A + F(u(t), t)
```

where ``u \\in \\R^n`` and  ``0 < \\varepsilon < 1``
``A`` must be a **periodic matrix** i.e. ``e^{t A} = e^{(t+\\pi) A}`` for any ``t \\in \\R``

## Argument :
- `prob::ODEProblem` : The problem to solve
- `alg::HOODEAB{order, ntau}` : the Adams-Bashforth HOODE algorithm

## Keywords :
- `dt` : duration of a time interval
- `kwargs...` : other keywords

"""
function SciMLBase.solve(
    prob::ODEProblem,
    alg::HOODEAB{order,ntau};
    dt = nothing,
    kwargs...,
) where {order,ntau}
    isSplitODEProblem(prob.problem_type) || error("HOODEAB alg need SplitODEProblem type")
    (!isnothing(dt) && haskey(kwargs, :nb_t)) &&
        error("Only one of dt and nb_t must be defined")
    haskey(kwargs, :order) &&
        error("order must defined as parameter of algorithm : HOODEAB(order)")
    linop = LinearHOODEOperator(prob.f.f1)
    fct = prob.f.f2.f
    p = typeof(prob.p) == SciMLBase.NullParameters ? missing : prob.p
    ho_prob = if haskey(prob.kwargs, :B)
        HOODEProblem(fct, prob.u0, prob.tspan, p, linop.A, linop.epsilon, prob.kwargs[:B])
    else
        HOODEProblem(fct, prob.u0, prob.tspan, p, linop.A, linop.epsilon)
    end
    if !haskey(kwargs, :nb_t) && !isnothing(dt)
        kwargs = (kwargs..., nb_t = Int64(round((prob.tspan[2] - prob.tspan[1]) / dt)))
    end

    sol = SciMLBase.solve(
        ho_prob,
        HOODETwoScalesAB();
        nb_tau = ntau,
        order = order,
        kwargs...,
    )
    other = OtherHOODE(sol.par_u0, sol.p_coef, sol.absprec, sol.relprec)
    return SciMLBase.build_solution(
        prob, # prob::Union{AbstractODEProblem,AbstractDDEProblem},
        HOODETwoScalesAB(), # alg,
        sol.t, # t,
        sol.u, # u,
        dense = sol.dense,
        interp = sol.interp,
        retcode = sol.retcode,
        stats = other,
    )
end

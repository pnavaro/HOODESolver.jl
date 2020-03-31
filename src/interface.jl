# module HiOscDiffEq
#using Reexport
# @reexport using DiffEqBase
using DiffEqBase

include("twoscales_pure_ab.jl")

struct HiOscDEProblem{T} <:DiffEqBase.DEProblem
    """The HiOscDE is `du/dt = (1/epsilon)*A*u + f(u,p,t)`."""
    f::Function
    """The initial condition is `u(tspan[1]) = u0`."""
    u0::Vector{T}
    """The solution `u(t)` will be computed for `tspan[1] ≤ t ≤ tspan[2]`."""
    tspan::Tuple{T,T}
    """Constant parameters to be supplied as the second argument of `f`."""
    p
    """Periodic Matrix of the problem"""
    A::Matrix
    """epsilon of the problem"""
    epsilon::T
end

# struct HiOscInterpolation{T} <: AbstractDiffEqInterpolation
#     t::T
#     u_caret::T
#     parphi::PreparePhi
#     order
# end
# function (interp::HiOscInterpolation)(t)
#     return _getresult(interp.u_caret, t, 
#     interp.parphi, 
#     interp.t[1], interp.t[end], 
#     interp.order)
# end
if VERSION >= v"1.3.1"
    abstract type AbstractHiOscSolution{T,N} <: DiffEqBase.AbstractTimeseriesSolution{T,N} end
else
    abstract type AbstractHiOscSolution{T,N} <: DiffEqBase.AbstractTimeseriesSolution{T,N,T} end
end
struct HiOscDESolution{T} <:AbstractHiOscSolution{T,T}
    u::Vector{Vector{T}}
    sol_u_caret::Union{Vector{Array{Complex{T},2}}, Missing}
    t::Vector{T}
    dense::Bool
    order::Integer
    parphi::PreparePhi
    prob::HiOscDEProblem{T}
    retcode
#    interp::HiOscInterpolation
    interp
    absprec
    relprec
end
function (sol::HiOscDESolution)(t)
    if sol.dense
        return _getresult(sol.sol_u_caret, t, sol.parphi, sol.t[1], sol.t[end], sol.order)
    else
        return undef
    end
end
# function DiffEqBase.build_solution{T}(prob::HiOscDEProblem{T}, 
#     sol::Vector{Vector{T}}, 
#     t::Vector{T}, 
#     fftsol::Vector{Array{T,2}}) where T<:AbstractFloat
#     return HiOscDESolution(sol, t, fftsol)
# end
# function DiffEqBase.build_solution{T}(prob::HiOscDEProblem{T}, 
#     sol::Vector{Vector{T}}, 
#     t::Vector{T}) where T<:AbstractFloat
#     return HiOscDESolution(sol, t, undef)
# end




function DiffEqBase.solve(prob::HiOscDEProblem{T}; 
    nb_tau=32, order=4, order_prep=order+2, dense=true, 
    nb_t=100, getprecision=dense
) where T<:AbstractFloat
    retcode = :Success
    if getprecision
        s1=DiffEqBase.solve(prob; 
    nb_tau=nb_tau, order=order, order_prep=order_prep, dense=dense, nb_t=nb_t, getprecision=false)
        s2=DiffEqBase.solve(prob; 
    nb_tau=nb_tau, order=order, order_prep=order_prep, dense=dense, nb_t=nb_t-1, getprecision=false)
        absprec = norm(s1.u[end] - s2.u[end])
        relprec = absprec/max( norm(s1.u[end]), norm(s2.u[end]))
        if dense
            for i = 1:10
                t = rand(typeof(prob.epsilon)) * (prob.tspan[2]-prob.tspan[1]) 
                    + prob.tspan[1]
                a, b = s1(t), s2(t)
                ap = norm(a-b)
                rp = ap/max(norm(a),norm(b))
                absprec = max(absprec,ap)
                relprec = max(relprec,rp)
            end
        end
        return HiOscDESolution(s1.u, s1.sol_u_caret, s1.t, dense, order,
                s1.parphi, prob, retcode, nothing, absprec*nb_t, relprec*nb_t)
    end 
    parphi = PreparePhi(prob.epsilon, nb_tau, prob.A, prob.f, t_0=prob.tspan[1], paramfct=prob.p)
    par_u0 = PrepareU0(parphi, order_prep, prob.u0)
    pargen = PrepareTwoScalePureAB(nb_t, prob.tspan[2], order, par_u0)
    sol, _, sol_u_caret = if dense
        twoscales_pure_ab(pargen, res_fft=true)
    else
        twoscales_pure_ab(pargen), undef, missing
    end
    return HiOscDESolution(
        reshape(mapslices(x->[x], sol, dims=1),size(sol,2)), 
        sol_u_caret, 
        collect(0:nb_t)*prob.tspan[2]/big(nb_t).+prob.tspan[1],
        dense,
        order,
        parphi, prob, retcode, nothing, undef, undef)
end


# end

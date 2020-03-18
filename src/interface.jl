module HiOscDiffEq
using Reexport
@reexport using DiffEqBase

include("twoscales_pure_ab.jl")

struct HiOscDEProblem{T} <:DiffEqBase.DEProblem
    """The HiOscDE is `du/dt = (1/epsilon)*A*u + f(u,p,t)`."""
    f::Function
    """The initial condition is `u(tspan[1]) = u0`."""
    u0::Vector{T}
    """The solution `u(t)` will be computed for `tspan[1] ≤ t ≤ tspan[2]`."""
    tspan::Vector{T}
    """Constant parameters to be supplied as the second argument of `f`."""
    p
    """Periodic Matrix of the problem"""
    A::Matrix
    """epsilon of the problem"""
    epsilon::T
end


struct HiOscDESolution{T} <:DiffEqBase.AbstractTimeseriesSolution
    sol::Vector{Vector{T}}
    sol_u_caret::Vector{Array{Complex{T},2}}
    t::Vector{T}
    dense::Bool
    order::Integer
    parphi::PreparePhi
    absprec
    relprec
end
function (sol::HiOscDESolution)(t)
    if sol.dense
        return _getresult(sol.sol_u_caret, t, sol.parphi, t[1], t[end], sol.order)
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
    nb_t=100, getprecision=true
) where T<:AbstractFloat
    if getprecision
        s1=DiffEqBase.solve(prob; 
    nb_tau=nb_tau, order=order, order_prep=order_prep, dense=dense, nb_t=nb_t, getprecision=false)
        s2=DiffEqBase.solve(prob; 
    nb_tau=nb_tau, order=order, order_prep=order_prep, dense=dense, nb_t=nb_t+1, getprecision=false)
        absprec = norm(s1.sol[end] - s2.sol[end])
        relprec = absprec/max( norm(s1.sol[end]), norm(s2.sol[end]))
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
        return HiOscDESolution(s1.sol, s1.sol_u_caret, s1.t, dense, order,
                s1.parphi, absprec, relprec)
    end 
    parphi = PreparePhi(prob.epsilon, nb_tau, prob.A, prob.f)
    par_u0 = PrepareU0(parphi, prob.u0, order_prep)
    pargen = PrepareTwoScalePureAB(nb_t, prob.tspan[2], order, par_u0;
    t_begin=prob.tspan[1])
    sol, sol_u_caret = if dense
         twoscales_pure_ab(pargen, res_fft=true)
    else
        twoscales_pure_ab(pargen), undef
    end
    sol = HiOscDESolution(sol, sol_u_caret, 
        collect(0:nb_t)*prob.tspan[2]/big(nb_t)+prob.tspan[1],
        dense,
        order,
        parphi, undef, undef)
end


end

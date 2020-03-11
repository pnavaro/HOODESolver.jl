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


struct HiOscDESolution{T} <:DiffEqBase.DESolution
    sol::Vector{Vector{T}}
    sol_u_chap::Vector{Array{Complex{T},2}}
    t::Vector{T}
    fftsol::Vector{Array{Complex{T},2}}
end
function DiffEqBase.build_solution{T}(prob::HiOscDEProblem{T}, 
    sol::Vector{Vector{T}}, 
    t::Vector{T}, 
    fftsol::Vector{Array{T,2}}) where T<:AbstractFloat
    return HiOscDESolution(sol, t, fftsol)
end
function DiffEqBase.build_solution{T}(prob::HiOscDEProblem{T}, 
    sol::Vector{Vector{T}}, 
    t::Vector{T}t) where T<:AbstractFloat
    return HiOscDESolution(sol, t, undef)
end

function DiffEqBase.solve(prob::HiOscDEProblem{T}; ) where T<:AbstractFloat


  return DiffEqBase.build_solution(sol, t)
end


end

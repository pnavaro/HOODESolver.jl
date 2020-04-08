using DiffEqBase


#using Reexport
# @reexport using DiffEqBase

include("twoscales_pure_ab.jl")

struct HiOscODEProblem{T} <:DiffEqBase.DEProblem
    """The HiOscODE is `du/dt = (1/epsilon)*A*u + f(u,p,t)`."""
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

struct HiOscInterpolation{T} <: DiffEqBase.AbstractDiffEqInterpolation
    t::Vector{T}
    u_caret::Vector{Array{Complex{T},2}}
    parphi::PreparePhi
    order
end
function (interp::HiOscInterpolation)(t)
    return _getresult(interp.u_caret, t, 
    interp.parphi, 
    interp.t[1], interp.t[end], 
    interp.order)
end

if typeof(DiffEqBase.AbstractTimeseriesSolution{Float64,Float64,Float64}) == DataType
    abstract type AbstractHiOscSolution{T,N} <: DiffEqBase.AbstractTimeseriesSolution{T,N,N} end
else
    abstract type AbstractHiOscSolution{T,N} <: DiffEqBase.AbstractTimeseriesSolution{T,N} end
end
struct HiOscODESolution{T} <:AbstractHiOscSolution{T,T}
    u::Vector{Vector{T}}
    t::Vector{T}
    dense::Bool
    order::Integer
    parphi::PreparePhi
    prob::HiOscODEProblem{T}
    retcode
    interp::Union{HiOscInterpolation, Nothing}
    absprec
    relprec
end
(sol::HiOscODESolution)(t) = sol.dense ? sol.interp(t) : undef
# function DiffEqBase.build_solution{T}(prob::HiOscODEProblem{T}, 
#     sol::Vector{Vector{T}}, 
#     t::Vector{T}, 
#     fftsol::Vector{Array{T,2}}) where T<:AbstractFloat
#     return HiOscODESolution(sol, t, fftsol)
# end
# function DiffEqBase.build_solution{T}(prob::HiOscODEProblem{T}, 
#     sol::Vector{Vector{T}}, 
#     t::Vector{T}) where T<:AbstractFloat
#     return HiOscODESolution(sol, t, undef)
# end




# function DiffEqBase.solve(prob::HiOscODEProblem{T};
"""
    function DiffEqBase.solve(prob::HiOscODEProblem{T}; 
    nb_tau::Integer=32, 
    order::Integer=4, 
    order_prep::Integer=order+2, 
    dense::Bool=true, 
    nb_t::Integer=100, 
    getprecision::Bool=dense) where T<:AbstractFloat

solver for Highly oscillatory problems, that an ODE of this form
``\\frac{\\delta u}{\\delta t} = \\frac{1}{\\varepsilon} A + F(u, t)``
where ``u \\in \\R^n`` and  ``0 < \\varepsilon < 1``

# Argument :
- `prob::HiOscODEProblem{T}` : The problem to solve

# Keywords :
- `nb_tau::Integer=32` : number of values of FFT transform, must be power of twoscales_pure_ab
- `order::Integer=4` : order of Adams-Bashforth method, and also of the interpolatation
- `order_prep::Integer=order+2` : order of the preparation
- `dense::Bool=true` : if true it is possible to compute solution at any time of the interval
- `nb_t::Integer=100` : number of period slices
- `getprecision::Bool=dense` : compute the absolute and relative precision

# Examples :
"""
function DiffEqBase.solve(prob::HiOscODEProblem{T}; 
    nb_tau::Integer=32, order::Integer=4, order_prep::Integer=order+2,dense::Bool=true, 
    nb_t::Integer=100, getprecision::Bool=dense
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
                t = rand(T) * (prob.tspan[2]-prob.tspan[1]) 
                    + prob.tspan[1]
                a, b = s1(t), s2(t)
                ap = norm(a-b)
                rp = ap/max(norm(a),norm(b))
                absprec = max(absprec,ap)
                relprec = max(relprec,rp)
            end
        end
        return HiOscODESolution(s1.u, s1.t, dense, order,
                s1.parphi, prob, retcode, s1.interp, absprec*nb_t, relprec*nb_t)
    end 
    parphi = PreparePhi(prob.epsilon, nb_tau, prob.A, prob.f, t_0=prob.tspan[1], paramfct=prob.p)
    par_u0 = PrepareU0(parphi, order_prep, prob.u0)
    pargen = PrepareTwoScalesPureAB(nb_t, prob.tspan[2], order, par_u0)
    u_mat, _, u_caret = if dense
        twoscales_pure_ab(pargen, res_fft=true)
    else
        twoscales_pure_ab(pargen), undef, undef
    end
    t = collect(0:nb_t)*(prob.tspan[2]-prob.tspan[1])/nb_t .+ prob.tspan[1]
    interp = dense ? HiOscInterpolation{T}(t, u_caret, parphi, order) : nothing
    return HiOscODESolution(
        reshape(mapslices(x->[x], u_mat, dims=1),size(u_mat,2)), 
        t,
        dense,
        order,
        parphi, prob, retcode, interp, undef, undef)
end

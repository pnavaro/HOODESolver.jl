using DiffEqBase


#using Reexport
# @reexport using DiffEqBase

include("twoscales_pure_ab.jl")

struct HiOscODEFunction{iip, N} <: DiffEqBase.AbstractODEFunction{iip}
    f
end
(fct::HiOscODEFunction{false,3})(u, p, t)=fct.f(u, p, t)
(fct::HiOscODEFunction{false,2})(u, p, t)=fct.f(u, p)
(fct::HiOscODEFunction{false,1})(u, p, t)=fct.f(u)
function (fct::HiOscODEFunction{true,4})(u, p, t)
    du = zero(u)
    fct.f(du, u, p, t)
    return du
end    


struct HiOscODEProblem{T} <:DiffEqBase.DEProblem
    """The HiOscODE is `du/dt = (1/epsilon)*A*u + f(u,p,t)`."""
    f::HiOscODEFunction
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
    """Matrix of linear problem to get the exact solution"""
    B::Union{Matrix,Missing}
    function HiOscODEProblem(
f, 
u0::Vector{T}, 
tspan::Tuple{T,T}, 
p, 
A::Matrix, 
epsilon::T,
B::Union{Matrix, Missing}
)  where {T}
        fct = if hasmethod(f, (Vector{T}, Vector{T}, Any, T) )
            HiOscODEFunction{true,4}(f)
        elseif hasmethod(f, (Vector{T}, Any, T) )
            HiOscODEFunction{false,3}(f)
        elseif hasmethod(f, (Vector{T}, Any) )
            HiOscODEFunction{false,2}(f)
        elseif hasmethod(f, (Vector{T},) )
            HiOscODEFunction{false,1}(f)
        else
            println("err !!!!!")
        end
        return new{T}(fct, u0, tspan, p, A, epsilon, B)
    end # end of function
end # end of struct
function HiOscODEProblem(f, u0, tspan, p, A, epsilon)
    return HiOscODEProblem(f, u0, tspan, p, A, epsilon, missing)
end
struct HiOscInterpolation{T} <: DiffEqBase.AbstractDiffEqInterpolation
    t::Vector{T}
    u_caret::Vector{Array{Complex{T},2}}
    parphi::PreparePhi
    order
end
function (interp::HiOscInterpolation)(t)
    return _getresult(interp.t, 
    interp.u_caret, t, 
    interp.parphi, 
    interp.t[1], interp.t[end], 
    interp.order)
#    return _getresult(interp.u_caret, t, 
#    interp.parphi, 
#    interp.t[1], interp.t[end], 
#    interp.order)
end
(interp::HiOscInterpolation)(vt::Vector{Float64})=interp.(vt)
(interp::HiOscInterpolation)(vt::Vector{BigFloat})=interp.(vt)
if typeof(DiffEqBase.AbstractTimeseriesSolution{Float64,Float64,Float64}) == DataType
    abstract type AbstractHiOscSolution{T,N} <: DiffEqBase.AbstractTimeseriesSolution{T,N,N} end
else
    abstract type AbstractHiOscSolution{T,N} <: DiffEqBase.AbstractTimeseriesSolution{T,N} end
end
struct HiOscODESolution{T} <:AbstractHiOscSolution{T,T}
    u::Vector{Vector{T}}
    u_tr::Union{Vector{Vector{T}}, Nothing}
    t::Vector{T}
    dense::Bool
    order::Integer
    par_u0::PrepareU0
    p_coef::CoefExpAB
    prob::HiOscODEProblem{T}
    retcode
    interp::Union{HiOscInterpolation, Nothing}
    tslocation
    absprec
    relprec
end
function HiOscODESolution(retcode::Symbol)
    return HiOscODESolution( undef, undef, undef, undef, undef, undef, 
    retcode, undef, undef, undef)
end
(sol::HiOscODESolution)(t) = sol.dense ? sol.interp(t) : undef

abstract type AbstractHiOscODEAlgorithm end

struct HiOscTwoScalesAB <: AbstractHiOscODEAlgorithm end
struct HiOscETDRK2 <: AbstractHiOscODEAlgorithm end
struct HiOscETDRK3 <: AbstractHiOscODEAlgorithm end
struct HiOscETDRK4 <: AbstractHiOscODEAlgorithm end


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
    getprecision::Bool=dense,
    verbose=100,
    par_u0::Union{PrepareU0,Missing}=missing,
    p_coef::Union{CoefExpAB,Missing}=missing
    ) where T<:AbstractFloat

solver for Highly oscillatory problems, that an ODE of this form
``\\frac{\\delta u(t)}{\\delta t} = \\frac{1}{\\varepsilon} A + F(u(t), t)``
where ``u \\in \\R^n`` and  ``0 < \\varepsilon < 1``
``A`` must be a periodic matrix i.e. ``e^{t A} = e^{(t+\\pi) A}`` for any ``t \\in \\R``

# Argument :
- `prob::HiOscODEProblem{T}` : The problem to solve

# Keywords :
- `nb_tau::Integer=32` : number of values of FFT transform, must be power of twoscales_pure_ab
- `order::Integer=4` : order of Adams-Bashforth method, and also of the interpolatation
- `order_prep::Integer=order+2` : order of the preparation
- `dense::Bool=true` : if true it is possible to compute solution at any time of the interval
- `nb_t::Integer=100` : number of period slices
- `getprecision::Bool=dense` : compute the absolute and relative precision
- `par_u0::Union{PrepareU0,Missing}=missing` : preparation data for u0
- `p_coef::Union{CoefExpAB,Missing}=missing` : coefficients for Adams-Bashforth method

# Examples :
"""
function DiffEqBase.solve(prob::HiOscODEProblem{T}; kwargs...) where T<:AbstractFloat
    return DiffEqBase.solve(prob, HiOscTwoScalesAB(); kwargs...)
end
function DiffEqBase.solve(prob::HiOscODEProblem{T}, alg::HiOscTwoScalesAB; 
    nb_tau::Integer=32, order::Integer=4, order_prep::Integer=order+2, dense::Bool=true, 
    nb_t::Integer=100, getprecision::Bool=dense, verbose=100,
    par_u0::Union{PrepareU0,Missing}=missing,
    p_coef::Union{CoefExpAB,Missing}=missing
) where T<:AbstractFloat
    retcode = :Success
    nb_tau = prevpow(2,nb_tau)
    traceln( 100,
    "solve function prob=$prob,\n nb_tau=$nb_tau, order=$order, order_prep=$order_prep, dense=$dense,\n nb_t=$nb_t, getprecision=$getprecision, verbose=$verbose";
    verbose=verbose,
)
    if getprecision
        s1=DiffEqBase.solve(prob; 
    nb_tau=nb_tau, order=order, order_prep=order_prep, dense=dense, nb_t=nb_t, getprecision=false, verbose=verbose-10, par_u0=par_u0)
        n_nb_t = Integer(floor(1.1373*nb_t+1))
        s2=DiffEqBase.solve(prob;
    nb_tau=nb_tau, order=order, order_prep=order_prep, dense=dense, nb_t=n_nb_t, getprecision=false, verbose=verbose-10, par_u0=s1.par_u0)
        absprec = norm(s1.u[end] - s2.u[end])
        relprec = absprec/max( norm(s1.u[end]), norm(s2.u[end]))
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
        return HiOscODESolution(s1.u, s1.u_tr, s1.t, dense, order,
                s1.par_u0, s1.p_coef, prob, retcode, s1.interp, 0,
                Float64(absprec), Float64(relprec))
    end 

    par_u0 = if ismissing(par_u0)
        parphi = PreparePhi(prob.epsilon, nb_tau, prob.A, prob.f, prob.B;
t_0=prob.tspan[1], paramfct=prob.p)
        PrepareU0(parphi, order_prep, prob.u0)
    else
        par_u0
    end
    pargen = PrepareTwoScalesPureAB(nb_t, prob.tspan[2], order, par_u0, 
p_coef=p_coef, verbose=verbose)
    if dense
        u_mat, _, u_caret = twoscales_pure_ab(pargen, res_fft=true, verbose=verbose)
        t = collect(0:nb_t)*(prob.tspan[2]-prob.tspan[1])/nb_t .+ prob.tspan[1]
        interp = HiOscInterpolation{T}(t, u_caret, par_u0.parphi, order)
        u = reshape(mapslices(x->[x], u_mat, dims=1),size(u_mat,2))
        u_tr = reshape(mapslices(x->[x], transpose(u_mat), dims=1),size(u_mat,1))
    else
        u = Array{Array{T,1},1}(undef,2)
        u[1] = copy(prob.u0)
        u[end] = twoscales_pure_ab(pargen, verbose=verbose, only_end=true)
        t = [prob.tspan[1], prob.tspan[2]]
        u_tr = nothing
        interp = nothing
    end
    return HiOscODESolution(
        u,
        u_tr,
        t,
        dense,
        order,
        par_u0, pargen.p_coef, prob, retcode, interp, 0, undef, undef)
end

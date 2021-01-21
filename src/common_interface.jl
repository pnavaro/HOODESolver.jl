using LinearAlgebra
using GenericSchur

include("interface.jl")
abstract type AbstractDiffEqOperator{T} end
abstract type AbstractDiffEqLinearOperator{T} <: AbstractDiffEqOperator{T} end

struct LinearHOODEOperator{T} <: AbstractDiffEqLinearOperator{T}
    epsilon::T
    A::Matrix
end
DiffEqBase.isinplace(linop::LinearHOODEOperator, n::Int)=false

function LinearHOODEOperator(mat::Matrix{T}) where{T}
    
    tab = zeros(Int,size(mat,1))
    sz = 0
    testlgn(x)=count(y->y!=0,x)!=0
    for i=1:size(mat,1)
        if testlgn(mat[i,:]) && testlgn(mat[:,i])
            sz += 1
            tab[sz] = i
        end
    end
    newmat = zeros(T,sz,sz)
    for i=1:sz, j=1:sz
        newmat[i,j] = mat[tab[i],tab[j]]
    end
    P=eigvecs(newmat)
    Ad = inv(P)*newmat*P

#    epsilon = 1/norm(Ad, Inf) ca bug avec Double64
    epsilon = 1/maximum(abs.(Ad))
    A = epsilon * mat
    LinearHOODEOperator{T}(epsilon, A)
end
function LinearHOODEOperator(linop::DiffEqArrayOperator)
    linop.update_func != DEFAULT_UPDATE_FUNC && error("no update operator function for HOODEAB Alg")
    LinearHOODEOperator(linop.A)
end
LinearHOODEOperator(linop::LinearHOODEOperator)=linop
LinearHOODEOperator(odefct::ODEFunction{false,LinearHOODEOperator{T}}) where{T}=odefct.f
isSplitODEProblem(probtype::Any)=false
isSplitODEProblem(probtype::SplitODEProblem)=true






struct HOODEAB{order, ntau} <: DiffEqBase.AbstractODEAlgorithm 
HOODEAB(order::Int=4; ntau=32)=new{order, ntau}()
end
export HOODEAB

function DiffEqBase.solve(prob::ODEProblem,
                          alg::HOODEAB{order, ntau};
                          dt = missing,
                          nb_t = missing) where{order, ntau}
    isSplitODEProblem(prob.problem_type) || error("HOODEAB alg need SplitODEProblem type")
    (!ismissing(dt) || !ismissing(nb_t)) && error("Only one of dt and nb_t must be defined")
    linop = LinearHOODEOperator(prob.f.f1)   
    fct = prob.f.f2.f
    p = typeof(prob.p) == DiffEqBase.NullParameters ? missing : prob.p
    ho_prob = HOODEProblem(fct, prob.u0, prob.tspan, p, linop.A, linop.epsilon)
    if ismissing(nb_t)
        nb_t = ismissing(dt) ? 100 : (Int)round((prob.tspan[2]-prob.tspan[1])/dt)
    end

    sol = DiffEqBase.solve(
        ho_prob,
        HOODETwoScalesAB();
        nb_tau = ntau,
        order = order,
        nb_t = nb_t
    )
    return sol
end



    


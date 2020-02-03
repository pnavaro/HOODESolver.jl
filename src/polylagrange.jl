#=
polylagrange:
- Julia version: 
- Author: ymocquar
- Date: 2019-11-25
=#

include("polyexp.jl")
using Polynomials

function getpolylagrange(k::Int64, j::Int64, N::DataType)
    @assert k <= j "_getpolylagrange(k=$k,j=$j) k must be less or equal to j"
    @assert N<:Signed "the type $N must be an Integer"
    result = Poly([one(Complex{Rational{N}})])
    for l=0:j
        if l != k
            result *= Poly([l//1,1//1])/(l-k)
        end
    end
    return result
end
struct CoefExpAB
    tab_coef
    function CoefExpAB(order::Int64, eps::AbstractFloat, lTau, dt::AbstractFloat)
        nTau = size(lTau,1)
        T = typeof(eps)
        tab_coef = zeros(Complex{T}, order+1, order+1, nTau)
        N = T == BigFloat ? BigInt : Int64
        eps = rationalize(N,eps, tol=eps*10*Base.eps(T) )
        dt = rationalize(N,dt, tol=dt*10*Base.eps(T) )
        lTau = rationalize.(N,lTau)
        pol_x = Poly([0//1,1//dt])
        for j=0:order
            for k=0:j
                res = view(tab_coef,k+1,j+1,:)
                pol = getpolylagrange(k, j, N)
                pol2 = pol(pol_x)
                for ind=1:nTau
                    ell = lTau[ind]
                    polint = if ell == 0
                        # in this case the exponentiel value is always 1
                        Polynomials.polyint(pol2)
                    else
                        polyint(PolyExp(pol2, im*ell/eps, -im*ell*dt/eps))
                    end
                    res[ind] = polint(dt)-polint(0)
                    if j==order && k in [4,5,6] && ind==1
                        println("j=$j k=$k ind=$ind pol=$pol")
                        println("j=$j k=$k ind=$ind pol2=$pol2")
                        println("j=$j k=$k ind=$ind polint=$polint")
                        println("2 j=$j k=$k ind=$ind res=$(res[ind])")
                    end
                end
            end
        end
        return new(tab_coef)
    end
end
struct CoefExpABfloat
    tab_coef
    function CoefExpABfloat(order::Int64, eps::AbstractFloat, lTau, dt)
        nTau = size(lTau,1)
        T=typeof(eps)
        N = T == BigFloat ? BigInt : Int64
        tab_coef = zeros(Complex{T}, order+1, order+1, nTau)
        pol_x = Poly([0, 1/dt])
        for j=0:order
            for k=0:j
                res = view(tab_coef,k+1,j+1,:)
                pol = getpolylagrange(k, j, N)
                pol2 = pol(pol_x)
                for ind=1:nTau
                    ell = lTau[ind]
                    polint = if ell == 0
                        # in this case the exponentiel value is always 1
                        Polynomials.polyint(pol2)
                    else
                        polyint(PolyExp(pol2, im*ell/eps, -im*ell*dt/eps))
                    end
                    res[ind] = polint(dt)-polint(0)
                    if j==order && k in [4,5,6] && ind==1
                        println("f0 j=$j k=$k ind=$ind pol=$pol")
                        println("f1 j=$j k=$k ind=$ind pol2=$pol2")
                        println("float j=$j k=$k ind=$ind polint=$polint")
                        println("fl2 j=$j k=$k ind=$ind res=$(res[ind])")
                    end


                end
            end
        end
        return new(tab_coef)
    end
end

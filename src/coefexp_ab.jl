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
    function CoefExpAB(order::Int64, epsilon::AbstractFloat, lTau, dt::AbstractFloat)
        nTau = size(lTau,1)
        T = typeof(epsilon)
        tab_coef = zeros(Complex{T}, order+1, order+1, nTau)
        N = T == BigFloat ? BigInt : Int64
        epsilon = rationalize(N,epsilon, tol=epsilon*10*Base.eps(T) )
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
                    pol3 = undef
                    pol_int = if ell == 0
                        # in this case the exponentiel value is always 1
                        Polynomials.polyint(pol2)
                    else
                        pol3=PolyExp(pol2, im*ell/epsilon, -im*ell*dt/epsilon)
                        polyint(pol3)
                    end
                    res[ind] = pol_int(dt)-pol_int(0)
 #                   if j==order && k in [3,5,7] && ind==2
                    #     println("epsilon=$epsilon dt=$dt")
                    #     println("j=$j k=$k ind=$ind pol=$pol")
                    #     println("j=$j k=$k ind=$ind pol2=$pol2")
                    #     println("j=$j k=$k ind=$ind pol3=$pol3")
                    #    println("j=$j k=$k ind=$ind pol_int=$pol_int")
 #                   println("pol_int(0)=$(pol_int(0))")
 #                   println("pol_int(dt)=$(pol_int(dt))")
                    #     println("2 j=$j k=$k ind=$ind res=$(res[ind])")
 #                   end

                end
            end
        end
        return new(tab_coef)
    end
end
struct CoefExpABfloat
    tab_coef
    function CoefExpABfloat(order::Int64, epsilon::AbstractFloat, lTau, dt)
        nTau = size(lTau,1)
        T=typeof(epsilon)
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
                    # res[ind] = if ell == 0
                    #     # in this case the exponentiel value is always 1
                    #     Polynomials.polyint(pol2)(dt)
                    # else
                    #     pol3 = polyint(PolyExp(pol2, im*ell/epsilon, -im*ell*dt/epsilon)).p
                    #     pol3(dt) - pol3.a[1]*exp(-im*ell*dt/epsilon)
                    # end
                    pol_int = if ell == 0
                        # in this case the exponentiel value is always 1
                        Polynomials.polyint(pol2)
                    else
                        polyint(PolyExp(pol2, im*ell/epsilon, -im*ell*dt/epsilon))
                    end
                    res[ind] = pol_int(dt)-pol_int(0)
#                    res[ind] = (ell == 0 ? pol_int(dt) : pol_int.p(dt)) - pol_int(0)
                    if j==order && k in [4,5,6] && ind==1
                        println("f0 j=$j k=$k ind=$ind pol=$pol")
                        println("f1 j=$j k=$k ind=$ind pol2=$pol2")
#                        println("float j=$j k=$k ind=$ind pol_int=$pol_int")
                        println("fl2 j=$j k=$k ind=$ind res=$(res[ind])")
                    end


                end
            end
        end
        return new(tab_coef)
    end
end


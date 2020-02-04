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
struct CoefExpABRational
    tab_coef
    function CoefExpABRational(order::Int64, epsilon::AbstractFloat, list_tau, dt::AbstractFloat)
        n_tau = size(list_tau,1)
        T = typeof(epsilon)
        tab_coef = zeros(Complex{T}, n_tau, order+1, order+1)
        N = T == BigFloat ? BigInt : Int64
        epsilon = rationalize(N,epsilon, tol=epsilon*10*Base.eps(T) )
        dt = rationalize(N,dt, tol=dt*10*Base.eps(T) )
        list_tau = rationalize.(N,list_tau)
        pol_x = Poly([0//1,1//dt])
        for j=0:order
            for k=0:j
                res = view(tab_coef, :, k+1, j+1)
                pol = getpolylagrange(k, j, N)
                pol2 = pol(pol_x)
                for ind=1:n_tau
                    ell = list_tau[ind]
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
struct CoefExpAB
    tab_coef
    function CoefExpAB(order::Int64, epsilon::AbstractFloat, n_tau, dt)
        T=typeof(epsilon)
        N, coef = T == BigFloat ? (BigInt, 10) : (Int64, 1)
        prec = 0
        eps_rat = rationalize(N, epsilon, coef*Base.eps(epsilon))
        dt_rat = rationalize(N, dt, coef*Base.eps(dt))
        if T == BigFloat
            prec = precision(BigFloat)
            setprecision(1024)
        end
        list_tau = [collect(0:n_tau / 2 - 1); collect(-n_tau / 2:-1)]
        epsilon = float(eps_rat)
        dt = float(dt_rat)
        T2 = BigFloat 
        tab_coef = zeros(Complex{T2}, n_tau, order+1, order+1)
        pol_x = Poly([0, 1/dt])
        for j=0:order
            for k=0:j
                res = view(tab_coef, :, k+1,j+1)
                pol = getpolylagrange(k, j, N)
                pol2 = pol(pol_x)
                for ind=1:n_tau
                    ell = list_tau[ind]
                    pol_int = if ell == 0
                        # in this case the exponentiel value is always 1
                        Polynomials.polyint(pol2)
                    else
                        polyint(PolyExp(pol2, im*ell/epsilon, -im*ell*dt/epsilon))
                    end
                    res[ind] = pol_int(dt)-pol_int(0)
                end
            end
        end
        if prec != 0
            setprecision(prec)
        end
        return new(T.(tab_coef))
    end
end


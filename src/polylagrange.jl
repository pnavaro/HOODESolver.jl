#=
PolyLagrangeOld:
- Julia version: 
- Author: ymocquar
- Date: 2019-11-25
=#
using Polynomials
import Base: show, *, +


function _getPolyLagrangeOld(k::Int64, j::Int64)
    @assert k <= j "_getPolyLagrangeOld(k=$k,j=$j) k must be less or equal to j"
    result = Poly([one(Complex{Rational{BigInt}})])
    for l=0:j
        if l != k
            result *= Poly([l//1,1//1])/(l-k)
        end
    end
    return result
end
function _getPolyLagrangeOld(k::Int64, j::Int64, dec::Int64)
    @assert k <= j && dec <= j "_getPolyLagrangeOld(k=$k,j=$j,dec=$dec) k and dec must be less or equal to j"
k_2 = (k+dec)%(j+1)-dec

    result = Poly([one(Complex{Rational{BigInt}})])
    for l=0:j
        if l != k
            l_2 = (l+dec)%(j+1)-dec
             result *= Poly([l_2//1,1//1])/(l_2-k_2)
        end
    end
    return result
end

# for j=0:10
#     for k=0:j
#         println("(k,j)=($k,$j) res=$(getPolyLagrangeOld(k,j))")
#     end
# end
struct PolyLagrangeOld
    poly
    polynum
#    polyexp
    respe
    respe_neg
    ressimple
    function PolyLagrangeOld(order::Int64, eps, lTau, dt)
      #  println("PolyLagrangeOld order=$order eps=$eps lTau[1:3]=$(lTau[1:3]), dt=$dt")
        zer = zero(Complex{Rational{BigInt}})
        poly = repeat([Poly([zer])], order+1, order+1)
        polynum = repeat([Poly([zer])], order+1, order+1, order+1)
        nTau = size(lTau,1)
        respe = zeros(Complex{BigFloat}, order+1, order+1, nTau)
        respe_neg = zeros(Complex{BigFloat}, order+1, order+1, nTau)
        ressimple = zeros(Rational{BigInt}, order+1, order+1, order+1)

         for j=0:order
            for k=0:j
                poly[k+1,j+1] = pol = _getPolyLagrangeOld(k, j)
                for z=0:j
                    polynum[k+1,j+1,z+1] = _getPolyLagrangeOld(k, j, z)
                end
            end
            for k=0:j
                res = view(respe,k+1,j+1,:)
                res_neg = view(respe_neg,k+1,j+1,:)
                pol = poly[k+1,j+1]
                pol2 = pol(Poly([0//1,1//dt]))
                pol2_neg = pol(Poly([0//1,-1//dt]))
                #               println("pol2=$pol2")
                for ind=1:nTau
                    ell = lTau[ind]
                    if ell == 0
                        # in this case the exponentiel value is always 1
                        polint = Polynomials.polyint(pol2)
                        polint_neg = Polynomials.polyint(pol2_neg)
 #                       println("k=$k, j=$j, ell=$ell, pol=$(polint) res=$(polint(dt)),$(convert(Complex{Float64},polint(dt)))")
 #                       println("k=$k, j=$j, ell=$ell, pol=$(convert(Poly{Complex{Float64}},polint))")
                   else
                    polint = polyint(PolyExp(pol2, big(im*ell/eps), big(-im*ell*dt/eps) ))
                    polint_neg = polyint(PolyExp(pol2_neg, big(im*ell/eps), big(im*ell*dt/eps) ))
                    if ind == 13 && (j == order || j == 6)
                        println("j=$j k=$k ind=$ind polint=$polint")
                        println("j=$j k=$k ind=$ind polint=$polint_neg")
                    end
                        #                       println("k=$k, j=$j, ell=$ell, pol=$(polint.p) res=$(polint.p(dt)),$(convert(Complex{Float64},polint.p(dt)))")
 #                      println("k=$k, j=$j, ell=$ell, pol=$(convert(Poly{Complex{Float64}},polint.p))")
                   end
                   res[ind] = polint(dt)-polint(0)
                   res_neg[ind] = polint_neg(-dt)-polint_neg(0)
                   if ind == 13 && (j == order || j == 6)
                    println("j=$j k=$k ind=$ind res=$(res[ind])")
                    println("j=$j k=$k ind=$ind res_neg=$(res_neg[ind])")
                   end

 #                  println("polint=$(string(polint))")
 #                  println("k=$k, j=$j, ell=$ell, res=$(convert(Complex{Float64},res[ind]))")
                end
                for z=0:j
                    p = polynum[k+1,j+1,z+1]
                    p_2 = p(Poly([0//1,1//dt]))
                    p_int = Polynomials.polyint(p_2)
                    ressimple[k+1,j+1,z+1] = (p_int(dt)-p_int(0))/dt
                end


 #                println("PolyLagrangeOld(k=$k,j=$j)=$(poly[k+1,j+1])")
            end
        end
        
  #      println("ressimple=$ressimple")

        return new(poly, polynum, respe, respe_neg, ressimple)
    end
end

# getPolyLagrangeOld( par::PolyLagrangeOld, k, j) = view(par.poly, k, j, 1:(j+1))
getPolyLagrangeOld( par::PolyLagrangeOld, k, j) = par.poly[k+1, j+1]
getPolyLagrangeOld( par::PolyLagrangeOld, k, j, dec) = par.polynum[k+1, j+1, dec+1]


"""
    PolyExp(p::Vector{T}, a::T, b::T)
    PolyExp(pol::Poly{T},a::T,b::T)

On the model of `Poly` from package `Polynomials`, construct a function that is a polynome multiply by an exponential
function. The exponential is an exponential of an affine function ``a x + b``.
The polynome is construct from its coefficients `p`, lowest order first.

If ``f = (p_n x^n + \\ldots + p_2 x^2 + p_1 x + p_0)e^{a x + b}``, we construct this through
`PolyExp([a_0, a_1, ..., a_n], a, b)`. 
It is also possible to construct it directly from the polynome.

In the sequels some methods with the same name than for Poly are implemented (`polyder`,
`polyint`, `strings`, ...) but not all, only the methods needed are developped.

# Arguments :
- `p::Vector{T}` or pol::Poly{T} : vector of coefficients of the polynome, or directly the polynome.
- `a::T`, `b::T` : coefficients of affine exponentiated function.

# Examples
```julia
julia> pe=PolyExp([1,2,3],2,1)
PolyExp(Poly(1 + 2*x + 3*x^2)*exp(2*x + 1))

julia> pe(0)
2.718281828459045

julia> pe(1)
120.51322153912601
```
"""
struct PolyExp{T}
    p::Poly{T}
    a::T
    b::T
    PolyExp(p::Vector{T},a::T,b::T) where{T<:Number}=new{T}(Poly{T}(p), a, b)
    PolyExp(pol::Poly{T},a::T,b::T) where{T<:Number}=PolyExp(Polynomials.coeffs(pol), a, b)
end
function _printNumberPar(x::Number) 
    return isreal(x) ? "$(real(x))" : (iszero(real(x)) ? "$(imag(x))im" : "($x)")
end
function _printNumber(x::Number)
    return isreal(x) ? "$(real(x))" : (iszero(real(x)) ? "$(imag(x))im" : "$x")
end
function Base.show(io::IO, pe::PolyExp)
    return print(
    io, 
    "PolyExp($(pe.p)*exp($(_printNumberPar(pe.a))*x + $(_printNumber(pe.b))))"
)
end
"""
    polyder(pe::PolyExp)

Construct the derivative of the `pe` function.

# Examples
```julia
julia> polyder(PolyExp([1, 3, -1],3,1))
PolyExp(Poly(6 + 7*x - 3*x^2)*exp(3*x + 1))

julia> polyder(PolyExp([1.0+im, 3im, -1, 4.0], 2.0+1.5im,1.0im))
PolyExp(Poly((0.5 + 6.5im) - (6.5 - 6.0im)*x + (10.0 - 1.5im)*x^2 + (8.0 + 6.0im)*x^3)*exp((2.0 + 1.5im)*x + 1.0im))
```
"""
function Polynomials.polyder(pe::PolyExp)
    return PolyExp(pe.p*pe.a + Polynomials.polyder(pe.p), pe.a, pe.b)
end

"""
    polyint(pe::PolyExp)

Construct the integrate function of `pe` which is of `PolyExp` type.
The algorithm used is a recursive integration by parts.

# Examples
```julia
julia> polyint(PolyExp([1.0,2,3],2.0,5.0))
PolyExp(Poly(0.75 - 0.5*x + 1.5*x^2)*exp(2.0*x + 5.0))

julia> polyint(PolyExp([1.0+0im,2],2.0im,3.0+0im))
PolyExp(Poly((0.5 - 0.5im) - 1.0im*x)*exp(2.0im*x + 3.0))
```
"""
function polyint(pe::PolyExp)
    if ( pe.a != 0 )
        pol = pe.p / pe.a
        if Polynomials.degree(pe.p) > 0
            pol -= polyint(PolyExp(Polynomials.polyder(pe.p), pe.a, pe.b)).p / pe.a
        end
    else
        pol = Polynomials.polyint(pol)
    end
    return PolyExp( pol, pe.a, pe.b)
end

polyval(pe::PolyExp, v::AbstractArray) = map(x->polyval(pe, x), v)

polyval(pe::PolyExp, x::Number) = pe.p(x) * exp(pe.a*x + pe.b)

(pe::PolyExp)(x) = polyval(pe, x)

coeffs(pe::PolyExp) = vcat(Polynomials.coeffs(pe.p), pe.a, pe.b)

function polyint(p::PolyExp, v_begin::Number, v_end::Number )
    pol = polyint(p)
    return pol(v_end)-pol(v_begin)
end

*( p1::PolyExp, p2::PolyExp )=PolyExp(p1.p*p2.p,p1.a+p2.a,p1.b+p2.b)

function +( p1::PolyExp, p2::PolyExp )
    @assert (p1.a == p2.a && p1.b == p2.b) "for adding exponents must be equal"
    return PolyExp(p1.p+p2.p,p1.a,p1.b)
end




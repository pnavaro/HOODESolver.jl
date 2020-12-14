using Polynomials
import Base: show, *, +



"""
    PolyExp(p::Vector{T}, a::T, b::T)
    PolyExp(pol::Polynomial{T},a::T,b::T)

On the model of `Polynomial` from package `Polynomials`, construct a function that is a polynome multiply by an exponential
function. The exponential is an exponential of an affine function ``a x + b``.
The polynome is construct from its coefficients `p`, lowest order first.

If ``f = (p_n x^n + \\ldots + p_2 x^2 + p_1 x + p_0)e^{a x + b}``, we construct this through
`PolyExp([p_0, p_1, ..., p_n], a, b)`. 
It is also possible to construct it directly from the polynome.

In the sequels some methods with the same name than for Polynomial are implemented
(`derivative`, `integrate`, `strings`, ...) but not all, only the methods needed are developped.

# Arguments :
- `p::Vector{T}` or pol::Polynomial{T} : vector of coefficients of the polynome, or directly the polynome.
- `a::T`, `b::T` : coefficients of affine exponentiated function.

# Examples
```julia
julia> pe=PolyExp([1,2,3],2,1)
PolyExp(Polynomial(1 + 2*x + 3*x^2)*exp(2*x + 1))

julia> pe(0)
2.718281828459045

julia> pe(1)
120.51322153912601
```
"""
struct PolyExp{T}
    p::Polynomial{T}
    a::T
    b::T
    PolyExp(p::Vector{T}, a::T, b::T) where {T<:Number} = new{T}(Polynomial{T}(p), a, b)
    PolyExp(pol::Polynomial{T}, a::T, b::T) where {T<:Number} = new{T}(pol, a, b)
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
        "PolyExp($(pe.p)*exp($(_printNumberPar(pe.a))*x + $(_printNumber(pe.b))))",
    )
end
"""
    derivative(pe::PolyExp)

Construct the derivative of the `pe` function.

# Examples
```julia
julia> derivative(PolyExp([1, 3, -1],3,1))
PolyExp(Polynomial(6 + 7*x - 3*x^2)*exp(3*x + 1))

julia> derivative(PolyExp([1.0+im, 3im, -1, 4.0], 2.0+1.5im,1.0im))
PolyExp(Polynomial((0.5 + 6.5im) - (6.5 - 6.0im)*x + (10.0 - 1.5im)*x^2 + (8.0 + 6.0im)*x^3)*exp((2.0 + 1.5im)*x + 1.0im))
```
"""
function Polynomials.derivative(pe::PolyExp)
    return PolyExp(pe.p * pe.a + Polynomials.derivative(pe.p), pe.a, pe.b)
end

"""
    integrate(pe::PolyExp)

Construct the integrate function of `pe` which is of `PolyExp` type.
The algorithm used is a recursive integration by parts.

# Examples
```julia
julia> integrate(PolyExp([1.0,2,3],2.0,5.0))
PolyExp(Polynomial(0.75 - 0.5*x + 1.5*x^2)*exp(2.0*x + 5.0))

julia> integrate(PolyExp([1.0+0im,2],2.0im,3.0+0im))
PolyExp(Polynomial((0.5 - 0.5im) - 1.0im*x)*exp(2.0im*x + 3.0))
```
"""
function Polynomials.integrate(pe::PolyExp)
    if (pe.a != 0)
        pol = pe.p / pe.a
        if Polynomials.degree(pe.p) > 0
            pol -=
                Polynomials.integrate(PolyExp(Polynomials.derivative(pe.p), pe.a, pe.b)).p /
                pe.a
        end
    else
        pol = Polynomials.integrate(pol)
    end
    return PolyExp(pol, pe.a, pe.b)
end
(pe::PolyExp)(x) = pe.p(x) * exp(pe.a * x + pe.b)
coeffs(pe::PolyExp) = vcat(Polynomials.coeffs(pe.p), pe.a, pe.b)
function Polynomials.integrate(p::PolyExp, v_begin::Number, v_end::Number)
    pol = Polynomials.integrate(p)
    return pol(v_end) - pol(v_begin)
end
*(p1::PolyExp, p2::PolyExp) = PolyExp(p1.p * p2.p, p1.a + p2.a, p1.b + p2.b)
function +(p1::PolyExp, p2::PolyExp)
    @assert (p1.a == p2.a && p1.b == p2.b) "for adding exponents must be equal"
    return PolyExp(p1.p + p2.p, p1.a, p1.b)
end

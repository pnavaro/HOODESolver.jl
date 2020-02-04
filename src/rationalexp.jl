import Base: show, *, convert

struct RationalExp
    coef::Complex{Rational{BigInt}}
    coef_exp::Complex{Rational{BigInt}}
end
*( re1::RationalExp, re2::RationalExp )=RationalExp(re1.coef*re2.coef,re1.coef_exp+re1.coef_exp)
function _printNumberPar(x::Number) 
    return isreal(x) ? "$(real(x))" : (iszero(real(x)) ? "$(imag(x))im" : "($x)")
end
function _printNumber(x::Number)
    return isreal(x) ? "$(real(x))" : (iszero(real(x)) ? "$(imag(x))im" : "$x")
end

function Base.show(io::IO, re::RationalExp)
    return print(
    io, 
    "RationalExp($(_printNumberPar(re.coef))*exp($(_printNumber(re.coef_exp)))"
)
end

Base.convert(BigFloat, x::RationalExp)=x.coef*exp(x.coef_par)


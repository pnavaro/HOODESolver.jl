using Base
using LinearAlgebra
using SparseArrays

function _expm2( mat )
    res = one(mat)
    resprec=zero(mat)
    mult = one(mat)
    i=1
    while resprec != res
        resprec .= res
        mult *= mat
        mult /= i
        res += mult
        i += 1
    end
 #   println("normdiff=$(norm(resprec-res))")
    return res
end


function _expm1( mat )
    valnorm= norm(mat)
    if valnorm > 1
        b = mat/valnorm
        coef_int = Integer(floor(valnorm))
        coef_mantissa = valnorm - coef_int
        return _expm2(b)^coef_int * _expm2(coef_mantissa*b)
    else
        return _expm2(mat)
    end
end

Base.exp(mat::Array{Complex{BigFloat}, 2}) = _expm1(mat)
Base.exp(mat::Array{BigFloat, 2}) = _expm1(mat)
Base.exp(mat::Array{Integer, 2}) = _expm1(mat)
Base.exp(mat::SparseMatrixCSC{Complex{BigFloat}, Integer}) = _expm1(mat)
Base.exp(mat::SparseMatrixCSC{BigFloat, Integer}) = _expm1(mat)
Base.exp(mat::Array{Rational, 2}) = Base.exp(float(mat))
Base.exp(mat::SparseMatrixCSC{Rational, Integer}) = Base.exp(float(mat))

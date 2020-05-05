using Base
using LinearAlgebra
using SparseArrays

function _expmat2( mat )
    res = one(mat)
    resprec = zero(mat)
    mult = one(mat)
    i=1
    cpt=0
    while cpt < 4
        resprec .= res
        mult *= mat
        mult /= i
        res += mult
        i += 1
        cpt = resprec == res ? cpt +1 : 0
    end
    return res
end
function _expmat1( mat )
    valnorm= norm(mat)

    if valnorm > 1
        b = mat/valnorm
        coef_int = Integer(floor(valnorm))
        coef_mantissa = valnorm - coef_int
        return _expmat2(b)^coef_int * _expmat2(coef_mantissa*b)
    else
        return _expmat2(mat)
    end
end
function _expmat0( mat )
    res  = setprecision(precision(BigFloat)+32) do
        _expmat1(mat)
    end
    return BigFloat.(res)
end
Base.exp(mat::Array{Complex{BigFloat}, 2}) = _expmat0(mat)
Base.exp(mat::Array{BigFloat, 2}) = _expmat0(mat)
Base.exp(mat::Array{Integer, 2}) = _expmat1(mat)
Base.exp(mat::SparseMatrixCSC{Complex{BigFloat}, Int64}) = _expmat0(mat)
Base.exp(mat::SparseMatrixCSC{BigFloat, Int64}) = _expmat0(mat)
Base.exp(mat::SparseMatrixCSC{Float64, Int64}) = _expmat1(mat)
Base.exp(mat::Array{Rational, 2}) = Base.exp(float(mat))
Base.exp(mat::SparseMatrixCSC{Rational, Integer}) = Base.exp(float(mat))

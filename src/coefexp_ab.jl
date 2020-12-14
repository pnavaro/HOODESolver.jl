#=
polylagrange:
- Julia version: 
- Author: ymocquar
- Date: 2019-11-25
=#
include("polyexp.jl")
function getpolylagrange(k::Int64, j::Int64, N::DataType)
    @assert k <= j "_getpolylagrange(k=$k,j=$j) k must be less or equal to j"
    @assert N <: Signed "the type $N must be an Integer"
    result = Polynomial([one(Complex{Rational{N}})])
    for l = 0:j
        if l != k
            result *= Polynomial([l // 1, 1 // 1]) / (l - k)
        end
    end
    return result
end
function interpolate(tab, order, value, N::DataType)
    T = (N == BigInt) ? BigFloat : Float64
    res = zeros(Complex{T}, size(tab[1]))
    for i = 0:order
        res .+= getpolylagrange(i, order, N)(-value) * tab[i+1]
    end
    return res
end
function getpolylagrange(t::Vector{T}, k, j) where {T<:Number}
    result = Polynomial([one(T)])
    for l = 0:j
        if l != k
            result *= Polynomial([-t[l+1], one(T)]) / (t[k+1] - t[l+1])
        end
    end
    return result
end
function interpolate(tab_time::Vector{T}, tab, order, value) where {T<:Number}
    res = zeros(Complex{T}, size(tab[1]))
    for i = 0:order
        res .+= getpolylagrange(tab_time, i, order)(value) * tab[i+1]
    end
    return res
end
interpolate(tab, order, value) = interpolate(tab, order, value, BigInt)
struct CoefExpABRational
    tab_coef::Any
    function CoefExpABRational(
        order::Int64,
        epsilon::AbstractFloat,
        list_tau,
        dt::AbstractFloat,
    )
        n_tau = size(list_tau, 1)
        T = typeof(epsilon)
        tab_coef = zeros(Complex{T}, n_tau, order + 1, order + 1)
        N = T == BigFloat ? BigInt : Int64
        epsilon = rationalize(N, epsilon, tol = epsilon * 10 * Base.eps(T))
        dt = rationalize(N, dt, tol = dt * 10 * Base.eps(T))
        list_tau = rationalize.(N, list_tau)
        pol_x = Polynomial([0 // 1, 1 // dt])
        for j = 0:order
            for k = 0:j
                res = view(tab_coef, :, k + 1, j + 1)
                pol = getpolylagrange(k, j, N)
                pol2 = pol(pol_x)
                for ind = 1:n_tau
                    ell = list_tau[ind]
                    pol3 = undef
                    pol_int = if ell == 0
                        # in this case the exponentiel value is always 1
                        Polynomials.integrate(pol2)
                    else
                        pol3 = PolyExp(pol2, im * ell / epsilon, -im * ell * dt / epsilon)
                        Polynomial.integrate(pol3)
                    end
                    res[ind] = pol_int(dt) - pol_int(0)
                end
            end
        end
        return new(tab_coef)
    end
end
struct CoefExpAB
    tab_coef::Any
    tab_coef_neg::Any
    function CoefExpAB(order::Int64, epsilon::AbstractFloat, n_tau, dt)
        T = typeof(epsilon)
        N = T == BigFloat ? BigInt : Int64
        new_prec = precision(BigFloat)
        new_prec += T == BigFloat ? order * 16 : 0
        setprecision(BigFloat, new_prec) do
            list_tau = [collect(0:n_tau/2-1); collect(-n_tau/2:-1)]
            T2 = BigFloat
            epsilon = T2(epsilon)
            dt = T2(dt)
            tab_coef = zeros(Complex{T2}, n_tau, order + 1, order + 1)
            pol_x = Polynomial([0, 1 / dt])
            for j = 0:order
                for k = 0:j
                    res = view(tab_coef, :, k + 1, j + 1)
                    pol = getpolylagrange(k, j, N)
                    pol2 = pol(pol_x)
                    for ind = 1:n_tau
                        ell = list_tau[ind]
                        pol_int = if ell == 0
                            # in this case the exponentiel value is always 1
                            Polynomials.integrate(pol2)
                        else
                            Polynomials.integrate(PolyExp(
                                pol2,
                                im * ell / epsilon,
                                -im * ell * dt / epsilon,
                            ))
                        end
                        res[ind] = pol_int(dt) - pol_int(0)
                    end
                end
            end # end of for j=....
        end # end of setprecision(...)
        # conversion to the new precision
        tab_coef = T.(real(tab_coef)) + im * T.(imag(tab_coef))
        tab_coef_neg = -conj(tab_coef)
        return new(tab_coef, tab_coef_neg)
    end
end

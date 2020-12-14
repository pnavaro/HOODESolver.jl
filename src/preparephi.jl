include("expmatrices.jl")
include("fftbig.jl")
using LinearAlgebra
using SparseArrays
using Statistics
"""
    PreparePhi(n_tau::Integer, epsilon::AbstractFloat, matrix_A::Matrix{Number},
    fct::Function, [matrix_B::Matrix])

Immutable structure, to share calculations, needed for the phi function.
These data can be used elsewhere for example in twoscale function.

# Arguments :
- `n_tau::Integer` : number of value around the unit disk, it must be a power of two.
- `epsilon::AbstractFloat` : epsilon of the system, the type of this value will be the typeof the result.
- `matrix_A::Matrix{Number}` : Matrix of twoscale system
- `fct::Function` : function of the system
- `[matrix_B::Matrix]` : matrix representing the linear function for debug

# Keywords
- `mode=1` : possibility of addionnal modes for optional behavior
- `paramfct=missing` : middle parameter of function fct
- `t_0=zero(epsilon)` : beginning of the time

# Fields :
- `epsilon` : epsilon of the system.
- `n_tau` : number of values for fft
- `tau_list` : list of values around the unit disk (0 ... n_tau/2 -n_tau/2-1 ... -1 )
- `tau_int` : list of values to integrate FFT
- `matrix_Ap` : sparse matrix with time dimension
- `tau_Ap` : for each angular ``\\tau`` around the unit disk the matrix ``e^{(\\tau \\time Ap)}``
- `tau_Ap_inv` : inverse of tau_Ap
- `par_fft` : fft parameters
- `fct` : function of differential equation
- `paramfct` : middle parameter of function fct
- `size_vect` : size of vector that is the size of the matrix
- `matrix_B` : B matrix for the linear case
- `mode` : for optional behavior
- `t_0` : beginning of the time

"""
struct PreparePhi
    epsilon::Any
    n_tau::Any
    tau_list::Any
    tau_int::Any
    sparse_Ap::Any
    tau_Ap::Any
    tau_Ap_inv::Any
    par_fft::Any
    par_fft_big::Any
    fct::Any
    paramfct::Any
    size_vect::Any
    matrix_B::Union{Matrix,Missing} # for debug only (linear function)
    mode::Any
    t_0::Any
    function PreparePhi(
        epsilon::AbstractFloat,
        n_tau::Integer,
        matrix_A::Matrix,
        fct::Function,
        matrix_B::Union{Matrix,Missing};
        mode = 1,
        paramfct = missing,
        t_0 = zero(epsilon),
    )
        T = typeof(epsilon)
        @assert prevpow(2, n_tau) == n_tau "$n_tau is not a power of 2"
        #        @assert isa(fct, Function) && hasmethod(fct, Tuple{Array{T}}) 
        #        "function fct is not correct"
        N = typeof(matrix_A[1, 1])
        size_vect = size(matrix_A, 1) + 1
        sparse_Ap = sparse(vcat(
            hcat(matrix_A, zeros(N, size_vect - 1)),
            transpose(zeros(N, size_vect)),
        ))
        tau_Ap = Vector{SparseMatrixCSC{T,Int64}}(undef, n_tau)
        tau_Ap_inv = Vector{SparseMatrixCSC{T,Int64}}(undef, n_tau)
        prec = precision(T)
        for i = 1:n_tau
            tau_Ap[i] =
                BigFloat.(
                    setprecision(prec + 32) do
                        round.(
                            exp((i - 1) * 2big(pi) / n_tau * sparse_Ap),
                            digits = prec + 16,
                            base = 2,
                        )
                    end,
                )
            tau_Ap_inv[i] =
                BigFloat.(
                    setprecision(prec + 32) do
                        round.(
                            exp(-(i - 1) * 2big(pi) / n_tau * sparse_Ap),
                            digits = prec + 16,
                            base = 2,
                        )
                    end,
                )
        end
        @assert tau_Ap[div(n_tau, 2)+1]^2 == I "The matrix must be periodic"
        tau_list = [collect(0:(div(n_tau, 2)-1)); collect(-div(n_tau, 2):-1)]
        tau_int = collect(transpose(vcat([0], -im ./ tau_list[2:end])))
        par_fft_big = PrepareFftBig(n_tau, big(epsilon))
        par_fft = T == BigFloat ? par_fft_big : missing
        return new(
            epsilon,
            n_tau,
            tau_list,
            tau_int,
            sparse_Ap,
            tau_Ap,
            tau_Ap_inv,
            par_fft,
            par_fft_big,
            fct,
            paramfct,
            size_vect,
            matrix_B,
            mode,
            t_0,
        )
    end
    function PreparePhi(
        epsilon::AbstractFloat,
        n_tau::Integer,
        matrix_A::Matrix,
        fct::Function;
        mode = 1,
        paramfct = missing,
        t_0 = zero(epsilon),
    )

        return PreparePhi(
            epsilon,
            n_tau,
            matrix_A,
            fct,
            missing,
            mode = mode,
            paramfct = paramfct,
            t_0 = t_0,
        )
    end
end

# filtred_f(u, mat_inv, fct, mat, p, t)= mat_inv * fct(mat*u, p, t)
# function filtredfct(par::PreparePhi, u_mat::Array{T,2}, t::T) where T <: Number
#  #   println("size=$(size(u_mat)) paramfct=$(par.paramfct)")
#     return reshape(
#     collect(Iterators.flatten(filtred_f.(eachcol(u_mat), par.tau_A_inv, par.fct, par.tau_A, par.paramfct, t))),
#     par.size_vect,
#     par.n_tau
# )
function homogeneousfct(par::PreparePhi, u)
    return vcat(par.fct(u[1:(end-1)], par.paramfct, u[end]), [1])
end
filtred_f(u, mat_inv, mat, par) = mat_inv * homogeneousfct(par, mat * u)
function filtredfct(par::PreparePhi, u_mat::Array{T,2}) where {T<:Number}
    #   println("size=$(size(u_mat)) paramfct=$(par.paramfct)")
    return reshape(
        collect(Iterators.flatten(filtred_f.(
            eachcol(u_mat),
            par.tau_Ap_inv,
            par.tau_Ap,
            (par,),
        ))),
        par.size_vect,
        par.n_tau,
    )
    # return filtred_f.(eachcol(u_mat), par.tau_A_inv, par.fct, par.tau_A)
end
isexactsol(par::PreparePhi) = !ismissing(par.matrix_B)
# for this function the time begins at par.t_begin
function getexactsol(par::PreparePhi, u0, t)
    @assert !ismissing(par.matrix_B) "The debug matrix is not defined"
    sparse_A = par.sparse_Ap[1:(end-1), 1:(end-1)]
    m = (1 / par.epsilon) * sparse_A + par.matrix_B
    t0 = par.t_0
    if ismissing(par.paramfct)
        return exp((t - t0) * m) * u0
    end
    a, b = par.paramfct
    mm1 = m^(-1)
    mm2 = mm1^2
    e_t0 = exp(-t0 * m)
    C = e_t0 * u0 + mm1 * e_t0 * (t0 * a + b) + mm2 * e_t0 * a
    e_inv = exp(-t * m)
    e = exp(t * m)
    C_t = -mm1 * e_inv * (t * a + b) - mm2 * e_inv * a
    return e * C + e * C_t
end

function phi(par::PreparePhi, u, order)
    @assert 2 <= order <= 20 "the order is $order, it must be between 2 and 20"
    f = undef
    if order == 2
        f = filtredfct(par, reshape(repeat(u, par.n_tau), par.size_vect, par.n_tau))
    else
        coef = if par.mode == 2
            eps(BigFloat)^0.5 # experimental
        else
            big(par.epsilon)^(order - 2)
        end
        #        coef = par.epsilon^(order/1.9117569711) #just to try
        #        coef = eps(typeof(par.epsilon))^0.2
        resPhi_u = phi(par, u, order - 1)
        f = resPhi_u + reshape(repeat(u, par.n_tau), par.size_vect, par.n_tau)
        f = filtredfct(par, f)
        #       f11 = coef * real(fftGen(par.par_fft, f)[:, 1]) / par.n_tau
        f11 = coef * mean(f, dims = 2)
        f .-= (phi(par, u + f11, order - 1) - resPhi_u) / coef
        if par.mode == 5
            f .-=
                par.epsilon^2 * reshape(
                    collect(Iterators.flatten((par.tau_A .- (I,)) .* (par.paramfct[1],))),
                    par.size_vect,
                    par.n_tau,
                )
        end

    end
    f = fftgen(par.par_fft_big, f)
    f = big(par.epsilon) * real(ifftgen(par.par_fft_big, f .* par.tau_int))
    return f
end
function get_tab_rand(T::DataType, s, n)
    tab = 2rand(T, s, n) - ones(T, s, n)
    correct = sum(tab, dims = 2) / n
    tab .-= correct
    return tab
end
"""
    PrepareU0(parphi::PreparePhi, order, u0, newprec)

Preparation of the original data

# Arguments
- `parphi::PreparePhi` : phi prepared parameters
- `order` : order of the preparation
- `u0` : initial data
- `[newprec]` : precision for the compute, if no given a default value is computed as a function of epsilon

# Fields
- `parphi` : phi prepared parameters
- `order` : order of the preparation
- `ut0` : formated initial data
- `u0` :initial data

"""
struct PrepareU0
    parphi::PreparePhi
    order::Any # it is the order of preparation at less one more than the order of AB process
    ut0::Any  # formated initial data
    up0::Any # initial data
    function PrepareU0(parphi::PreparePhi, order, u0, newprec)
        #
        # Numerical preparation (corresponds to formula (2.6), p3 with j=2, and
        # algorithm is given p4).
        up0 = vcat(u0, [parphi.t_0])
        up0 = big.(up0)
        #    y, um =
        #     if parphi.mode == 3
        #         tab = transpose(reshape(repeat(get_tab_rand(parphi.n_tau, 
        # typeof(parphi.epsilon)), parphi.size_vect), parphi.n_tau, parphi.size_vect))
        #         y = parphi.epsilon*tab
        #         um = up0 - y[:,1]
        #         y, um
        #     else
        if newprec == 0
            prec = precision(BigFloat)
            newprec = prec + 32 + div(-exponent(parphi.epsilon) * order^2, 3)
            #            newprec = prec*4
            #            println("prec : $prec --> $newprec")
        end
        y, um = setprecision(newprec) do
            y = phi(parphi, up0, 2)
            um = up0 - y[:, 1]
            for i = 3:order
                y = phi(parphi, um, i)
                um = up0 - y[:, 1]
            end
            y, um
        end
        #      end
        ut0 = reshape(repeat(um, parphi.n_tau), parphi.size_vect, parphi.n_tau) + y
        # if parphi.mode == 4
        #     ut0 +=parphi.epsilon^2*reshape(collect(Iterators.flatten((parphi.tau_A .- (I,)).* (parphi.paramfct[1],))),parphi.size_vect,parphi.n_tau)
        # end
        if typeof(parphi.epsilon) != BigFloat
            ut0 = convert.(Float64, ut0)
            up0 = convert.(Float64, up0)
        end
        return new(parphi, order, ut0, up0)
    end
    PrepareU0(parphi::PreparePhi, order, u0) = PrepareU0(parphi, order, u0, 0)
end

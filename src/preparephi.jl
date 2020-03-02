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

# Implementation :
- epsilon : epsilon of the system.
- n_tau : number of values for fft
- list_tau : list of value around the unit disk
- matrix_A : sparse matrix
- tau_A : for each angular ``\tau`` around the unit disk the matrix ``e^{(\tau \time A)}```
- par_fft : fft parameters.
- fct : function of differential equation.
- size_vect : size of vector that is the size of the matrix

"""
struct PreparePhi
    epsilon
    n_tau
    tau_list
    tau_int
    sparse_A 
    tau_A
    tau_A_inv
    par_fft
    fct
    size_vect
    mode
    matrix_B::Union{Matrix,Missing} # for debug only (linear function)
    function  PreparePhi(
    epsilon::AbstractFloat, 
    n_tau::Integer, 
    matrix_A::Matrix, 
    fct::Function,
    matrix_B::Union{Matrix,Missing};
    mode=1
)
        T = typeof(epsilon)
        @assert prevpow(2,n_tau) == n_tau "$n_tau is not a power of 2"
        @assert isa(fct, Function) && hasmethod(fct, Tuple{Array{T}}) 
        "function fct is not correct"
        sparse_A = sparse(matrix_A)
        tau_A = Vector{ SparseMatrixCSC{T,Int64}}(undef, n_tau)
        tau_A_inv = Vector{ SparseMatrixCSC{T,Int64}}(undef, n_tau)
        prec= precision(T)
        for i = 1:n_tau
            tau_A[i] = BigFloat.( setprecision(prec+32) do
                round.( 
    exp((i-1) * 2big(pi) / n_tau * sparse_A), 
    digits=prec+16, 
    base=2 
)
            end )
            tau_A_inv[i] = BigFloat.( setprecision(prec+32) do
                round.( 
    exp(-(i-1) * 2big(pi) / n_tau * sparse_A), 
    digits=prec+16, 
    base=2 
)
            end )
        end
        @assert tau_A[div(n_tau, 2) + 1]^2 == I "The matrix must be periodic"
        size_vect = size(matrix_A,1)
        tau_list = [collect(0:(div(n_tau, 2) - 1)); collect(-div(n_tau, 2):-1)]
#         tau_list_2 = collect(transpose(reshape(repeat(
#     tau_list,
#     size_vect
# ),n_tau,size_vect)))       
        tau_int = collect(transpose(vcat([0], -im ./ tau_list[2:end])))
        par_fft = T == BigFloat ? PrepareFftBig(n_tau, epsilon) : missing
        return new( 
    epsilon, 
    n_tau, 
    tau_list,
    tau_int,
    sparse_A, 
    tau_A, 
    tau_A_inv, 
    par_fft,
    fct,
    size_vect,
    matrix_B,
    mode
)
    end
    function  PreparePhi(
        epsilon::AbstractFloat, 
        n_tau::Integer, 
        matrix_A::Matrix, 
        fct::Function;
        mode=1
    )
        return PreparePhi(epsilon, n_tau, matrix_A,fct, missing, mode)
    end
end
isexactsol(par::PreparePhi) = !ismissing(par.matrix_B)
function getexactsol(par::PreparePhi, u0, t)
    @assert !ismissing(par.matrix_B) "The debug matrix is not defined"
    return exp(t*(1/par.epsilon*par.sparse_A+par.matrix_B))*u0
end

filtred_f(u, mat_inv, fct, mat)= mat_inv * fct(mat*u)
function filtredfct(par::PreparePhi, u_mat::Array{T,2}) where T <: Number
    return reshape(
    collect(Iterators.flatten(filtred_f.(eachcol(u_mat), par.tau_A_inv, par.fct, par.tau_A))),
    par.size_vect,
    par.n_tau
)
# return filtred_f.(eachcol(u_mat), par.tau_A_inv, par.fct, par.tau_A)
end
function phi( par::PreparePhi, u, order)
    @assert 2 <= order <= 20 "the order is $order, it must be between 2 and 20"
    f = undef
    if  order == 2
        f = filtredfct(par, reshape(repeat(u, par.n_tau), par.size_vect, par.n_tau))
    else

        coef = if par.mode == 1
            par.epsilon^(order - 2)
        elseif par.mode == 2
            eps(BigFloat)^0.5
        end
#        coef = par.epsilon^(order/1.9117569711) #just to try
#        coef = eps(typeof(par.epsilon))^0.2
        resPhi_u = phi(par, u, order - 1)
        f = resPhi_u + reshape(repeat(u, par.n_tau), par.size_vect, par.n_tau)
        f = filtredfct(par, f)
 #       f11 = coef * real(fftGen(par.par_fft, f)[:, 1]) / par.n_tau
        f11 = coef * mean(f, dims=2)
        f .-= (phi(par, u + f11, order - 1) - resPhi_u) / coef
    end
    f = fftgen(par.par_fft, f)
    f = par.epsilon*real(ifftgen(par.par_fft, f .* par.tau_int))
    return f
end
struct PrepareU0
    parphi
    order # it is the order of preparation at less one more than the order of AB process
    ut0  # formated initial data
    u0 # initial data
    function PrepareU0(parphi::PreparePhi, order, u0, newprec)
        #
        # Numerical preparation (corresponds to formula (2.6), p3 with j=2, and
        # algorithm is given p4).
        if newprec == 0
            prec = precision(BigFloat)
            newprec = prec + 32 + div(-exponent(parphi.epsilon)*order^2, 3)
#            newprec = prec*4
            println("prec : $prec --> $newprec")
        end
        y, um = setprecision(newprec) do
            y = phi(parphi, u0, 2)
            um = u0 - y[:, 1]
            for i=3:order
                y = phi(parphi, um, i)
                um = u0 - y[:, 1]
            end
            y, um
        end
        ut0 = reshape(repeat(um, parphi.n_tau), parphi.size_vect, parphi.n_tau) + y
        return new(parphi, order, ut0, u0)
    end
    PrepareU0(parphi::PreparePhi, order, u0)=PrepareU0(parphi, order, u0, 0)
end

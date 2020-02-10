include("expmatrices.jl")
include("fftbig.jl")
using LinearAlgebra
using SparseArrays
using Statistics
"""
    PreparePhi(n_au::Integer, epsilon::AbstractFloat, matrix_A::Matrix{Number},
    fct::Function)
)

Immutable structure, to share calculations, needed for the phi function.
These data can be used elsewhere for example in twoscale function.

# Arguments :
- `n_tau::Integer` : number of value around the unit disk, it must be a power of two.
- `epsilon::AbstractFloat` : epsilon of the system, the type of this value will be the typeof the result.
- `matrix_A::Matrix{Number}` : Matrix of twoscale system
- `fct::Function` : function of the system

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
    epsilon3, 
    n_tau, 
    tau_list, 
    sparse_A, 
    tau_A, 
    par_fft,
    fct2,
    size_vect
    function  PreparePhi(
    n_tau::Integer, 
    epsilon4::AbstractFloat, 
    matrix_A::Matrix{Number}, 
    fct1::Function
)
        T = typeof(epsilon)
        @assert prevpow(2,n_tau) == n_tau "$n_tau is not a power of 2"
        @assert isa(fct, Function) && hasmethod(fct, Tuple{Array{T}, T}) 
        "function fct is not correct"
        sparse_A = sparse(matrix_A)
        tau_A = Vector{ SparseMatrixCSC{T,Int64}}(undef, n_tau)
        prec= precision(T)
        for i = 1:n_tau
            tau_A[i] = BigFloat.( setprecision(prec+32) do
                round.( 
    exp((i-1) * 2big(pi) / n_tau * sparse_A), 
    digits=prec+16, 
    base=2 
)
            end )
        end
        @assert tau_A[div(n_tau, 2) + 1]^2 == I "The matrix must be periodic"
        tau_list = [collect(0:(div(n_tau, 2) - 1)); collect(-div(n_tau, 2):-1)]       
        tau_int = vcat( [0], -im / tau_list[2,end]) # Coefficients to integrate
        par_fft = T == BigFloat ? PrepareFftBig(n_tau, epsilon) : missing
        return new( 
    epsilon, 
    n_tau, 
    tau_list, 
    sparse_A, 
    tau_A, 
    par_fft,
    fct,
    size(matrix_A,1),
)
    end
end
function filtredfct(par::PreparePhi, u_mat::Array{T,2}) where T <: Number
    return reshape(
    collect(Iterators.flatten(par.fct.( par.tau_A .* eachcol(u_mat), par.epsilon))),
    par.size_vect,
    par.n_tau
)
end
function phi( par::PreparePhi, u, order)
    @assert 2 <= order <= 20 "the order is $order, it must be between 2 and 20"
    if  order == 2
        f = filtredfct(par, reshape(repeat(u, par.n_tau), par.size_vect, par.n_tau))
    else
        coef = par.epsilon^(order - 2)
        resPhi_u = phi(par, u, order - 1)
        f = resPhi_u + reshape(repeat(u, par.n_tau), par.size_vect, par.n_tau)
        f = filtredfct(par, f)
 #       f11 = coef * real(fftGen(par.par_fft, f)[:, 1]) / par.n_tau
        f11 = coef * mean(f, dims=2)
        f .-= (phi(par, u + f11, order - 1) - resPhi_u) / coef
    end
    f = fftgen(par.parFft, f)
    f = par.epsilon*real(ifftgen(par.parFft, f .* par.tau_int))
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
        prec = 0
        if newprec != 0
            prec = precision(BigFloat)
            setprecision(newprec)
        end
        y = phi(2, u0, parphi)
        um = u0 - y[:, 1]
        for i=3:order
            y = phi(i, um, parphi)
            um = u0 - y[:, 1]
        end
        if prec != 0
            setprecision(prec)
        end
        ut0 = reshape(repeat(um, parphi.n_tau), parphi.size_vect, parphi.n_tau) + y
        return new(parphi, order, ut0, u0)
    end
    PrepareU0(parphi::PreparePhi, order, u0)=PrepareU0(parphi, order, u0, 0)
end

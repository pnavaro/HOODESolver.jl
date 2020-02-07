include("expmatrices.jl")
include("fftbig.jl")
using LinearAlgebra
using SparseArrays

"""
    PreparePhi(nTau::Integer, eps::AbstractFloat)

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
    epsilon, 
    n_tau, 
    tau_list, 
    sparse_A, 
    tau_A, 
    par_fft
    fct::Function
    size_vect
    function  PreparePhi(
    n_tau::Integer, 
    epsilon::AbstractFloat, 
    matrix_A::Matrix{Number}, 
    fct::Function
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
function filtredFctGen(par::PreparePhi, u_mat::Array{T,2}) where T <: Number
    return reshape(
    collect(Iterators.flatten(par.fct.( par.tau_A .* eachcol(u_mat), par.epsilon))),
    par.size_vect,
    par.n_tau
)
end

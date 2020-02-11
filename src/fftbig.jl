#=
fft:
- Julia version: 
- Author: ymocquar
- Date: 2019-11-15
=#
# function that reverse the order of the pos lowest bits
using FFTW
function _reverse_num(num, pos)
    result = 0
    pos_m1 = pos-1
    for i=0:pos_m1
        if (num & (1 << i))  != 0
           result |= 1 << (pos_m1 - i)
        end
    end
    return result
end
"""
    PrepareFftBig( size::Unsigned, [T=BigFloat])

Immutable structure to operate fft transform, 
x is the type of non transformed data also called signal.

# Arguments :
- `size::Integer` : Number of values, must be a power of two
- `[T=BigFloat | x::T ]` : type of the values

# Implementation
- size : size of the signal
- tab_permut : permutation
- root_one : size order roots of one
- root_one_conj : conjugate of root_one

"""
struct PrepareFftBig
    size
    tab_permut
    root_one
    root_one_conj
    type::DataType
    function PrepareFftBig( size::Integer, x::T ) where {T<:AbstractFloat}
        @assert prevpow(2,size) == size "size=$size is not a power of 2"
        power = convert(Int64,log2(size))
        tab_permut = zeros(Int64,size)
        root_one = zeros(Complex{T}, size )
        root_one_conj = zeros(Complex{T}, size )
        for i=1:size
            tab_permut[i] = _reverse_num( i-1, power) + 1
            root_one[i] = exp(one(T)*2pi*im*(i-1)/size)
            root_one_conj[i] = exp(-one(T)*2pi*im*(i-1)/size)
        end
        return new(size, tab_permut, root_one, root_one_conj, typeof(x))
    end
end
PrepareFftBig( size::Integer, type::DataType ) = PrepareFftBig( size, one(type))
PrepareFftBig( size::Integer) = PrepareFftBig( size, one(BigFloat))
function fftbig!(par::PrepareFftBig, signal; flag_inv=false)
    s=size(signal,2)
    @assert prevpow(2,s) == s "size(signal)=$s is not a power of 2"
    s_div2 = div(s,2)
    len = s
    n_len = len>>1
    nb_r = 1
    rootO = flag_inv ? par.root_one : par.root_one_conj;
    while n_len != 0
        start = 1
        suite = start+n_len
          for i=1:nb_r
 #           println("start=$start suite=$suite n_len_m1=$n_len_m1")
             deb = 1
             for j=deb:nb_r:s_div2
                signal[:,start], signal[:,suite] = (signal[:,start] + signal[:,suite,]), 
    (signal[:,start] - signal[:,suite])*rootO[j]
 #               println("start=$start suite=$suite nb_r=$nb_r")
                start += 1
                suite += 1
             end
             start = suite
             suite = start+n_len
        end
        len = n_len
        n_len >>= 1
        nb_r <<= 1
    end
    signal .= signal[:,par.tab_permut]       
    
    if flag_inv
        signal ./= s
    end
    return signal
end
function fftbig(par::PrepareFftBig, signal; flag_inv=false)
    fl = flag_inv
    return fftbig!(
        par::PrepareFftBig,
        copy(convert(Array{Complex{par.type}}, signal));
        flag_inv=fl
    )
end
fftgen(_::Any, t::Array{Complex{Float64}}) = fft(t, (2,))
fftgen(_::Any, t::Array{Float64}) = fft(t, (2,))
fftgen(p::PrepareFftBig, t::Array{T}) where {T<:AbstractFloat} = fftbig(p, t)
fftgen(p::PrepareFftBig, t::Array{Complex{T}}) where {T<:AbstractFloat} = fftbig(p, t)
ifftgen(_::Any, t::Array{Complex{Float64}}) = ifft(t, (2,))
ifftgen(_::Any, t::Array{Float64}) = ifft(t, (2,))
ifftgen(p::PrepareFftBig, t::Array{T}) where {T<:AbstractFloat}  = fftbig(p, t, flag_inv = true)
ifftgen(p::PrepareFftBig, t::Array{Complex{T}}) where {T<:AbstractFloat}  = fftbig(p, t, flag_inv = true)

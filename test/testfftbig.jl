#=
testfftbig:
- Julia version: 
- Author: ymocquar
- Date: 2019-11-18
=#
using Test
using FFTW
using Random
using LinearAlgebra

include("../src/fftbig.jl")

function getfalse( tab )
    for i=1:size(tab,1), j=1:size(tab,2)
        if !tab[i,j]
            return i,j
        end
    end
    return 0, 0
end



function testfftbig( s, T::DataType, seed_val )

    Random.seed!(seed_val)
    tab = zeros(Complex{T}, 1, s)
    tab .= rand(T, 1, s)
    if T == Float64
        tabfftref = fft(tab,2)
    else
        tab2 = zeros(Complex{Float64}, 1, s)
        tab2 .= tab
        tabfftref = fft(tab2,2)
    end

    tab_test = copy(tab)

    p = PrepareFftBig(s, real(tab[1, 1]))

    fftbig!(p, tab_test)

    @test isapprox(tabfftref, tab_test, atol=1e-15, rtol=1e-15)

    tol = (T == BigFloat) ? 1e-50 : 1e-15

    fftbig!(p, tab_test, flag_inv=true)

    @test getfalse(isapprox.(tab, tab_test, atol=tol, rtol=tol)) == (0, 0)

    @test isapprox(tab, tab_test, atol=tol, rtol=tol)

end
function testfftbig2( s, T::DataType, seed_val, nb_v )

    Random.seed!(seed_val)
    tab = zeros(Complex{T}, nb_v, s)
    tab .= rand(T, nb_v, s)
    if T == Float64
        tabfftref = fft(tab,2)
    else
        tab2 = zeros(Complex{Float64}, nb_v, s)
        tab2 .= tab
        tabfftref = fft(tab2,2)
    end

    tab_test = copy(tab)

    p = PrepareFftBig(s, one(T) )

    tab_test2 = fftbig(p, tab_test)

    @test isapprox(tabfftref, tab_test2, atol=1e-15, rtol=1e-15)

    fftbig!(p, tab_test)

    @test isapprox(tabfftref, tab_test, atol=1e-15, rtol=1e-15)

    tol = (T == BigFloat) ? 1e-50 : 1e-15

    tab_test3 = fftbig(p, tab_test, flag_inv=true)
    @test isapprox(tab, tab_test3, atol=tol, rtol=tol)

    fftbig!(p, tab_test, flag_inv=true)

    @test isapprox(tab, tab_test, atol=tol, rtol=tol)

end

function testfftbigprec(s, seed_v)

    @time @testset "test fftbig changing precision s=$s seed=$seed_v" begin
        Random.seed!(seed_v)
        for k=1:10
            setprecision(k*256) do
                p = PrepareFftBig(s)
                tab = rand(BigFloat, 5, s)
                tabRef = Complex.(tab)
                tabRes = fftgen(p,tab)
                tabRes2 = ifftgen(p, tabRes)
                println("k=$k norm=$(norm(tabRef-tabRes2))")
#                @test isapprox(tabRef, tabRes2, atol=1e-78^k)
#                @test isapprox(real.(tabRef), real.(tabRes2), atol=1e-78^k)
                @test isapprox(tabRef, tabRes2, atol=1e-75^k)
                @test isapprox(real.(tabRef), real.(tabRes2), atol=1e-75^k)
            end
        end
    end

end


tab_decl = [[8, Float64, 12345678], [8, BigFloat, 9876],[2^8, Float64, 1928], [2^10, BigFloat, 5656]]

for t in tab_decl
    s = t[1]
    type = t[2]
    seed_v = t[3]
    @time @testset "test fftbig for value size=$s type=$type" begin testfftbig(s, type, seed_v) end
end
tab_decl2 = [[8, Float64, 4556789, 4], [8, BigFloat, 4563, 4],[2^10, Float64, 9900771, 4], [2^13, BigFloat, 125609, 4]]

for t in tab_decl2
    s = t[1]
    type = t[2]
    seed_v = t[3]
    nb_v = t[4]
    @time @testset "test fftbig for vector size=$s type=$type nb_v=$nb_v" begin testfftbig2(s, type, seed_v, nb_v) end
end

tab_decl3 =[ 8965, 1919191, 56188827, 9017555]
for sd in tab_decl3
    testfftbigprec(32,sd)
end

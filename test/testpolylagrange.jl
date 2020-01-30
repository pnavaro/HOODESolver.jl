#=
testPolyLagrange:
- Julia version: 
- Author: ymocquar
- Date: 2019-11-25
=#
using Test


include("../src/polylagrange.jl")
#include("../src/twoScale3.jl")
#include("../out/dataTest_4.jl")
#include("dataTest_4.jl")

# function fct( x, poly)
#     val = 1.0
#     res = 0.0
#     for coef in poly
#         res += val*coef
#         val *= x
#     end
#     return res
# end


function testPolyLagrange()
    orderMax=10
    
        ind = 1
        eps = 1//1000
        nTau = 32
        nMax = 100
        tMax = 1//1

        @time @testset "test polylagrange order max=$orderMax eps=$eps nTau=$nTau" begin
             par = PolyLagrange(orderMax, eps, [collect(0:nTau ÷ 2 - 1); collect(-nTau ÷ 2:-1)], tMax//nMax)

            polyOri = Poly(convert( Array{Complex{Rational{BigInt}}},
                                    [3//4, 13//24, 1//7, 2//9, 40//51,
#                                     1//128, 2//79, 3//334, 7//19, 1//100,
#                                     3//5, 13//23, 1//8, 1//9, 41//51,
                                     1//144, 3//32, 7//19, 2//5, 1//100]
                                    ) )
            s = size(polyOri,1)
            polyResult = Poly([zero(Complex{Rational{BigInt}})])

            for i = 1:s
        #        x = convert(Float64,1-i)
                x = 1-i
                polyResult += polyOri(x)*getpolylagrange(par,i-1,s-1)
            end
      #     println("polyResult=$polyResult")
      #   @test isapprox(polyOri, polyResult, atol=1e-12, rtol=1e-12)
            @test polyOri == polyResult

            for dec = 0:(s-1)
     #           println("dec=$dec")
                polyResult = Poly([zero(Complex{Rational{BigInt}})])
                for i = 1:s
                    x =  - ((i-1+dec)%s-dec)
                    polyResult += polyOri(x)*getPolyLagrange(par,i-1,s-1, dec)
                end
                @test polyOri == polyResult
            end



    #         d = PreparePhi( nTau, Float64(eps) )


    #         p_ts = PrepareTwoScale(nMax, Float64(tMax), d)

    # #         println("pl=$(p_ts.pl)")
    # #         println("ql=$(p_ts.ql)")
    # #         println("p1=$(p_ts.p1)")
    # #         println("p2=$(p_ts.p2)")
    # #         println("p3=$(p_ts.p3)")
    #         prec=1e-12
    #         @test isapprox(p_ts.p1, par.respe[1,3,:], atol=prec, rtol=prec)
    #         @test isapprox(p_ts.p2, par.respe[2,3,:], atol=prec, rtol=prec)
    #         @test isapprox(p_ts.p3, par.respe[3,3,:], atol=prec, rtol=prec)

    #         pnl = tab[ind]; ind+= 1
    #         qnl = tab[ind]; ind+= 1
    #         p_0 = tab[ind]; ind+= 1
    #         p_1 = tab[ind]; ind+= 1
    #         p_2 = tab[ind]; ind+= 1
    #         p_3 = tab[ind]; ind+= 1
    #         prec=1e-10
    #         @test isapprox(p_0, par.respe[1,4,:], atol=prec, rtol=prec)
    #         @test isapprox(p_1, par.respe[2,4,:], atol=prec, rtol=prec)
    #         @test isapprox(p_2, par.respe[3,4,:], atol=prec, rtol=prec)
    #         @test isapprox(p_3, par.respe[4,4,:], atol=prec, rtol=prec)

    end
end

function testPolyExp()

    @time @testset "test PolyExp" begin
        pol = PolyExp(convert( Vector{Complex{Rational{BigInt}}}, [2//3, 5//23, 57//4, 1//5]), big(im//7), big(1//10+0im) )

        c_res = coeffs(pol)

 #       println("pol=$(string(pol))")

        @test [2//3, 5//23, 57//4, 1//5, im//7, 1//10 ] == c_res

        polder = polyder(pol)

        c_resder = coeffs(polder)

        @test [2im//21+5//23, 5im//(23*7)+ 57//2, 57im//(4*7)+3//5, 1im//35,  im//7, 1//10 ] == c_resder

        polintder = polyint(polder)
        c_resintder = coeffs(polintder)

        @test c_res == c_resintder

        p1 = PolyExp([1.0, 2, 3], 3.4, 1.1)
        p2 = PolyExp([0.5, -1, 1], 1.3, 0.245)

        c_resmul = coeffs(p1*p2)
        @test isapprox([0.5, 0.0, 0.5, -1.0, 3.0, 4.7, 1.345], c_resmul, atol=1e-15)

        p3 = PolyExp([0.2, 0, 3.0, 2.0], 1.3, 0.245)

        c_resadd = coeffs(p2 + p3)
        @test isapprox([0.7, -1.0, 4.0, 2.0, 1.3, 0.245], c_resadd, atol=1e-15)
    end
end

testPolyLagrange()

testPolyExp()

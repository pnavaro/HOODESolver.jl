
include("../src/expmatrices.jl")
using LinearAlgebra
using Random
using Test
function testexp()
    A=zeros(Int64, 4, 4)
    A[1,3] = 1
    A[3,1] = -1
    @time @testset "test 1 exp big matrices" begin
        @test isapprox( one(A), _expm1(0A), atol=1e-16)
        @test isapprox(exp(0.123*A), _expm1(0.123*A), atol=1e-15)
    end    
    @time @testset "test 2 exp big matrices" begin
        Abig = convert(Array{BigFloat,2},A)
        Random.seed!(81725)
        sens=1
        for i = 1:100:100000
            v=rand(BigFloat)*i*sens
            res = _expm0(v*Abig)
            resRef = one(Abig)
            resRef[1,1] = resRef[3,3] = cos(v)
            resRef[1,3] = sin(v)
            resRef[3,1]= -resRef[1,3]
 #           println("i=$i norm=$(norm(resRef-res))")
            @test isapprox( resRef, res, atol=1e-78)
            sens *= -1
        end
    end    
    @time @testset "test 3 exp big matrices" begin
        b = [ 0 0 0 1; 0 0 -1 0; 0 1 0 0; -1 0 0 0]
        for i=1:20
            C=rand(4,4)*i
            B = b*i
            @test isapprox(exp(C), _expm1(C), rtol=1e-13)
            @test isapprox(exp(21C), _expm1(21C), rtol=1e-12)
            @test isapprox(exp(22C), _expm1(22C), rtol=1e-12)
            @test isapprox(exp(23C), _expm1(23C), rtol=1e-11)
            @test isapprox(exp(24C), _expm1(24C), rtol=1e-10)
            @test isapprox(exp(B), _expm1(B), rtol=1e-14)
            @test isapprox(exp(10B), _expm1(10B), rtol=1e-13)
            @test isapprox(exp(100B), _expm1(100B), rtol=1e-12)
            @test isapprox(exp(1000B), _expm1(1000B), rtol=1e-11)
            @test isapprox(exp(10000B), _expm1(10000B), rtol=1e-10)
        end
    end    
    @time @testset "test 4 exp big matrices big precision" begin
        prec_old = prec = precision(BigFloat)
        for i=1:6
            prec *= 2
            setprecision(prec)
            @test isapprox(_expm1(2big(pi)*A)-I,zeros(BigFloat,4,4),atol=eps(BigFloat)*100)
        end
        setprecision(prec_old)
    end
end
testexp()

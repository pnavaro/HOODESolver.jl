
include("../src/expmatrices.jl")
using LinearAlgebra
using Random
using Test
function testexp()
    A=zeros(Int64, 4, 4)
    A[1,3] = 1
    A[3,1] = -1
    @time @testset "test 1 exp big matrices" begin
        @test isapprox( one(A), _expmat1(0A), atol=1e-16)
        @test isapprox(exp(0.123*A), _expmat1(0.123*A), atol=1e-15)
    end    
    @time @testset "test 2 exp big matrices" begin
        Abig = convert(Array{BigFloat,2},A)
        Random.seed!(81725)
        sens=1
        for k=1:4
            setprecision(k*256) do
                for i = 1:100:10000
                    v=rand(BigFloat)*i*sens
                    res = _expmat0(v*Abig)
                    resRef = one(Abig)
                    resRef[1,1] = resRef[3,3] = cos(v)
                    resRef[1,3] = sin(v)
                    resRef[3,1]= -resRef[1,3]
        #           println("i=$i norm=$(norm(resRef-res))")
                    @test isapprox( resRef, res, atol=(1e-77)^k)
                    sens *= -1
                end
            end
        end
    end    
    @time @testset "test 3 exp big matrices" begin
        Random.seed!(8781115)

        b = [ 0 0 0 1; 0 0 -1 0; 0 1 0 0; -1 0 0 0]
        for i=1:20
            C=rand(4,4)*i
            B = b*i
            @test isapprox(exp(C), _expmat1(C), rtol=1e-13)
            @test isapprox(exp(0.7C), _expmat1(0.7C), rtol=1e-12)
            @test isapprox(exp(0.89C), _expmat1(0.89C), rtol=1e-12)
            @test isapprox(exp(B), _expmat1(B), rtol=1e-14)
            @test isapprox(exp(10B), _expmat1(10B), rtol=1e-13)
            @test isapprox(exp(100B), _expmat1(100B), rtol=1e-12)
            @test isapprox(exp(1000B), _expmat1(1000B), rtol=1e-11)
            @test isapprox(exp(10000B), _expmat1(10000B), rtol=1e-10)
        end
    end    
    @time @testset "test 4 exp big matrices big precision" begin
        Random.seed!(90908)
        prec_old = prec = precision(BigFloat)
        for i=1:6
            prec *= 2
            setprecision(prec)
            @test isapprox(_expmat1(2big(pi)*A)-I,zeros(BigFloat,4,4),atol=eps(BigFloat)*100)
        end
        setprecision(prec_old)
    end
end
testexp()

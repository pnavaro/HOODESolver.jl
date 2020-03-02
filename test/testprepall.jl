include("../src/preparephi.jl")
include("../src/henon_heiles.jl")
using Test

function testprepare()
    n_tau=32
    A = [0 0 1 0; 0 0 0 0;-1 0 0 0; 0 0 0 0]
    u =rand(BigFloat, 4)
    for i=1:5
        prec = 256*2^(i-1)
        setprecision(256*2^(i-1))
        for i_eps = 1:40
            epsilon = big"1.25"^(-i_eps)
            println("prec=$prec epsilon=$epsilon")
            for order=3:10
                @time @testset "test PrepareU0 prec=$prec order=$order epsilon=$epsilon" begin
                    newprec = prec + 32 + div(-exponent(epsilon)*order^2, 3)
                    precref = newprec + 10
                    parphi = PreparePhi(epsilon, n_tau, A, henon_heiles)
                    pref=PrepareU0(parphi, order, u, precref)
                    pnew=PrepareU0(parphi, order, u, newprec)
                    @test pref.ut0 == pnew.ut0
                end
            end
        end
    end
end
testprepare()

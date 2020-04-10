include("../src/interface.jl")
using LinearAlgebra
using Test

function testcoeffloat()
    tab=[ (big"0.1", big"0.0001"), (big"0.00001", big"0.01")]
    @testset "test coefexp float vs bigfloat" begin
        for (bigeps, bigdt) in tab
            epsilon = convert(Float64,bigeps)
            dt = convert(Float64, bigeps)
            bigceab = CoefExpAB(10,bigeps,32,bigdt)
            ceab = CoefExpAB(10,epsilon,32,bigdt)
            @test isapprox(bigceab.tab_coef, ceab.tab_coef, atol=1e-15, rtol=1e-14)
    #       println("epsilon=$epsilon norm=$(norm(bigceab.tab_coef-ceab.tab_coef, Inf))")
        end
    end
 end
function testcomp()
    @testset "test solve float vs bigfloat" begin

        bigeps = big"0.0001"
        epsilon = convert(Float64,bigeps)
        bigu0=rand(BigFloat,4)
        u0 = convert.(Float64, bigu0)
        A =  [0 0 1 0; 0 0 0 0;-1 0 0 0; 0 0 0 0]
        fct = (u,p,t) -> [0, u[4], -2u[1]*u[2], -u[2]-u[1]^2+u[2]^2] # Henon-Heiles
        bigprob = HiOscODEProblem(fct, bigu0, (big"0.0",big"3.0"), missing, A, bigeps)
        prob = HiOscODEProblem(fct, u0, (0.0, 3.0), missing, A, epsilon)
        bigsol = solve(bigprob)
        sol = solve(prob) 

        for i=1:101 
            @test isapprox(
    bigsol.interp.u_caret[i],
    sol.interp.u_caret[i], 
    atol=1e-13, 
    rtol=1e-12
)
        end
    end
end

testcoeffloat()
testcomp()

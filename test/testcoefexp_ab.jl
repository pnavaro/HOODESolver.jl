include("../src/coefexp_ab.jl")
include("datacoefexp.jl")
using LinearAlgebra
using Test
function testcoefexp_ab() 
    par = CoefExpAB(15, big"0.1", 32, big"0.0001")
    tab, list_j = get_coef_ab_for_test()
    @time @testset "test coef for exponential Adams-Bashforth" begin
        for i in list_j
            println("norm = $(norm(par.tab_coef[:,:,j]-tab[:,:,j], Inf))")
            @time @test all(isapprox.(par.tab_coef[:,:,j], tab[:,:,j], atol=1e-76, rtol=1e-76))
        end
    end
end
testcoefexp_ab()

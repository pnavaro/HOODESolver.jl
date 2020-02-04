include("../src/coefexp_ab.jl")
include("datacoefexp.jl")
using LinearAlgebra
using Test
function testcoefexp_ab()
    
    par = CoefExpAB(15, big"0.1", 32, big"0.0001")

    tab = get_coef_ab_for_test()

    println("norm = $(norm(par.tab_coef-tab, Inf))")

    @test isapprox.(par.tab_coef, tab, atol=1e-70, rtol=1e-70)
end

testcoefexp_ab()

include("../src/coefexp_ab.jl")
using Test


function testcoefexp_ab()
    
    par = CoefExpAB(15, big"0.1", )

    @test 

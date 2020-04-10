include("../src/coefexp_ab.jl")
using LinearAlgebra

function testcoeffloat()
    tab=[ (big"0.1", big"0.0001"), (big"0.00001", big"0.01")]

    for (bigeps, bigdt) in tab
        epsilon = convert(Float64,bigeps)
        dt = convert(Float64, bigeps)
        bigceab = CoefExpAB(10,bigeps,32,bigdt)
        ceab = CoefExpAB(10,epsilon,32,bigdt)
        print("epsilon=$epsilon norm=$(norm(bigceab.tab_coef-ceab.tab_coef, Inf))")
    end
 end
 testcoeffloat()

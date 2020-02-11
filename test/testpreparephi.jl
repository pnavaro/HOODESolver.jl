include("../src/preparephi.jl")
include("../src/henon_heiles.jl")
include("dataprep_u0.jl")
using Test
function testpreparephi()
    @time @testset "phi and prepare u0" begin
        tab_ref = setprecision(1024) do
            get_prepare_u0_for_test()
        end
        for i_prec=1:3
            prec = 256*2^(i_prec-1)
            tol=max((1e-77)^i_prec, 1e-200)
            setprecision(prec) do
                epsilon = big"1"/256
                n_tau = 32
                parphi = PreparePhi(
    epsilon, 
    n_tau, 
    [0 0 1 0; 0 0 0 0; -1 0 0 0; 0 0 0 0],
    henon_heiles
)
                u0 = BigFloat.([0.125, 0.140625, 0.15625, 0.171875])
                for ord=2:10
            #     paru0 = PrepareU0(parphi, ord, u0, 1024)
                    paru0 = PrepareU0(parphi, ord, u0)
 #                   println("prec=$prec ord=$ord norm=$(norm(tab_ref[:, :, ord]- paru0.ut0))")
                    @test isapprox(tab_ref[:, :, ord], paru0.ut0, atol=tol, rtol=tol)
                end
            end
        end
    end
end
testpreparephi()

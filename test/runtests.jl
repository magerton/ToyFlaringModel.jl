using Revise
using ToyFlaringModel
using Optim
using Calculus
using Test

using ToyFlaringModel: G, Gp, F, Fp, A, Ap,
    cstar, cstarp, profit, k0, Va, cf,
    GatherCorners

@testset "flaring model derivatives etc" begin
    m = FlaringModel()
    mc = FlaringModel(monop=false)
    @test derivative(x -> G(m,x,0), 1.0) ≈ Gp(m, 1, 0)
    @test derivative(x -> F(m,x,0), 1.0) ≈ Fp(m, 1, 0)
    @test derivative(x -> A(m,x,0), 1.0) ≈ Ap(m, 1, 0)
    @test Gp(m,1,0) + Fp(m,1,0) + Ap(m,1,0) == 0

    @test Cost(mc,0) == k0(m)
    @test Va(mc,k0(mc),1) == Va(mc,k0(mc),0)

    @test Total(m,0) == Total(m,1) == Total(m,0) == Total(m,1)

    find_cstar(m,t) = optimize(x -> -profit(m,x,t), k0(m), 5)
    test_cstar(m,t) = Optim.minimizer(find_cstar(m,t))

    @test test_cstar(m,1) ≈ cstar(m,1)
    @test test_cstar(m,0) ≈ cstar(m,0)

    @test cstarp(m,0) ≈ derivative(x -> cstar(m,x), 0.0)
    @test cstarp(m,1) ≈ derivative(x -> cstar(m,x), 1.0)

    @test GatherCorners(m, k0(m), 0)[3][1] == Va(m, k0(m), 0)
    @test GatherCorners(m, k0(m), 0) isa Vector{NTuple{2,Float64}}
end


# cstarp(m,0)
#
#
# using ToyFlaringModel: a, d, b, e
#
# let t = 0,
#     x = cstar(m,t),
#     aa = a(m,t),
#     dd = d(m,t),
#     bb = b(m,t),
#     ee = e(m,t)
#
#     @show 3x^2 - 2(aa+bb)x + (aa*bb + dd*ee)
#     y = (2x - (aa+dd))/(6x - 2(aa+bb))
#     @show y, cstarp(m,t)
# end

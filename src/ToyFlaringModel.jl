module ToyFlaringModel

using Plots: Shape

"creates function f(x::T) = x.f for f in flds"
macro getFieldFunction(T, flds::Symbol...)
    r = quote end
    for f in flds
        e = esc( quote $f( x::$T ) = x.$f end )
        push!(r.args, e)
    end
    return r
end

export FlaringModel

struct FlaringModel
    Ro::Float64    # oil revenue
    Rg::Float64    # gas revenue
    cbar::Float64  # max flaring cost
    Vbar::Float64  # max value of alternate
    k0::Float64    # MC of processing (const)
    monop::Bool
end

# Access fields of flaring model
@getFieldFunction FlaringModel Ro Rg cbar Vbar k0 monop

# default constructor
FlaringModel(;Ro=2, Rg=1, cbar=4, Vbar=4, k0=1.25, monop=true) =
    FlaringModel(Ro,Rg,cbar,Vbar,k0,monop)

# total revenues
Rt(m) = Ro(m) + Rg(m)

# Thresholds
Va(m,c,t) = Rt(m) - c
cf(m,c,t) = c - Rg(m) - t

# Gather, Flare, AltOp and derivatives wrt c
G( m,c,t) = Va(m,c,t)*(cbar(m) - cf(m,c,t))
Gp(m,c,t) = -(Va(m,c,t) + cbar(m) - cf(m,c,t))
F( m,c,t) = cf(m,c,t)*(Va(m,c,t) + cf(m,c,t)/2)
Fp(m,c,t) = Va(m,c,t)
A( m,c,t) = cbar(m)*(Vbar(m)-Va(m,c,t)) - cf(m,c,t)^2/2
Ap(m,c,t) = cbar(m) - cf(m,c,t)

# Profit & profit prime
profit( m,c,t) = (c-k0(m))*G(m,c,t)
profitp(m,c,t) = (c-k0(m))*Gp(m,c,t) + G(m,c,t)

# some helper functions to compue cstar
a(m,t) = k0(m)
d(m,t) = Rt(m)
e(m,t) = cbar(m) + Rg(m) + t
b(m,t) = d(m,t) + e(m,t)

apb(m,t) = a(m,t)+b(m,t)
ab( m,t) = a(m,t)*b(m,t)
de( m,t) = d(m,t)*e(m,t)

cstar( m,t) = (apb(m,t) - sqrt( inside_sqrt(m,t) ) )/3
inside_sqrt(m,t) = apb(m,t)^2 - 3*(ab(m,t) + de(m,t))

# inside_sqrtp(m,t) = 2*apb(m,t) - 3*(a(m,t) + d(m,t))
# cstarp(m,t) = 1/3 - inside_sqrtp(m,t)/sqrt(inside_sqrt(m,t))/6

function cstarp(m,t)
    x = cstar(m,t)
    aa = a(m,t)
    dd = d(m,t)
    bb = b(m,t)
    dxdt = (2x - (aa+dd)) / (6x - 2(aa+bb))
    return dxdt
end


export Cost, Gather, Flare, AltOpt, Total

Cost(  m, t) = monop(m) ? cstar(m,t) : k0(m)
Gather(m, t) = G(m, Cost(m,t), t)
Flare( m, t) = F(m, Cost(m,t), t)
AltOpt(m, t) = A(m, Cost(m,t), t)
Total( m, t) = Gather(m,t) + Flare(m,t) + AltOpt(m,t)


function GatherCorners(m,cg,t)
    c = cf(m,cg,t)
    cb = cbar(m)
    V = Va(m,cg,t)
    o = 0.0
    return [ (o,c), (o,cb), (V, cb), (V, c) ]
end

function FlareCorners(m,cg,t)
    c = cf(m,cg,t)
    V = Va(m,cg,t)
    o = 0.0
    return [ (o, o), (o, c), (V, c), (V+c, o) ]
end

function AltOptCorners(m,cg,t)
    c = cf(m,cg,t)
    cb = cbar(m)
    V = Va(m,cg,t)
    Vb = Vbar(m)
    o = 0.0
    return [(V,cb), (Vb,cb), (Vb,o), (V+c,o), (V,c)]
end

export GatherShape, FlareShape, AltOptShape

GatherShape(m,c,t) = Shape( GatherCorners(m, c, t) )
FlareShape( m,c,t) = Shape( FlareCorners( m, c, t) )
AltOptShape(m,c,t) = Shape( AltOptCorners(m, c, t) )

GatherShape(m,t) = GatherShape( m, Cost(m,t), t)
FlareShape( m,t) = FlareShape(  m, Cost(m,t), t)
AltOptShape(m,t) = AltOptShape( m, Cost(m,t), t)

center(s::Shape) = sum(s.x)/length(s.x), sum(s.y)/length(s.y)

export GatherLab, FlareLab, AltOptLab

GatherLab(m,c,t) = (center(GatherShape(m,c,t))..., "Gather")
FlareLab( m,c,t) = (center(FlareShape( m,c,t))..., "Flare")
AltOptLab(m,c,t) = (center(AltOptShape(m,c,t))..., "Alternative")

GatherLab(m,t) = GatherLab(m, Cost(m,t), t)
FlareLab( m,t) = FlareLab( m, Cost(m,t), t)
AltOptLab(m,t) = AltOptLab(m, Cost(m,t), t)

end # module

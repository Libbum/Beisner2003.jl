using DifferentialEquations
using Plots

# Parameters
const f = 2.0 # depensation in piscivores (B)
const m = 0.01 # yr⁻¹ natural piscivore mortality
const a = 0.37 # B recruitment assymptote (shape parameter)
const r = 8.7 # mg⋅m⁻²⋅d⁻¹ maximum recycling of P
const q = 8 # Or 2, steepness coefficient of R
const k = 90 # mg⋅m⁻² half saturation constant of P recycling
const s = 0.001 # d⁻¹ flushing rate
const b = 0.5 # adult density when recruitment = a/2 (B recruitment shape parameter)
const α = 0.8 # mg⋅m⁻²⋅d⁻¹ maximum growth rate of A per P input rate
const v = 0.07 # m⋅d⁻¹ sinking rate per algal cell
const g = 0.1 # d⁻¹ C grazing coefficient

# Drivers
# H yr⁻¹ harvest rate
H = 0.1
# DOC mg/m2 dissolved organic carbon
DOC = 4.74 # Exp(1.556)
# L mg⋅m⁻²⋅d⁻¹ P loading
L = 1.429; # Exp(0.357)

# State variables
# B g/m2 piscivore population biomass (pike)
## W g/m2 planktivore population biomass (bream)
#W = @. 3.2*exp(-3.6926*B)
## C mm zooplankton community length
#C = W < 2.58 ? 0.5 : 0.35 # mm
## Z m depth of thermocline
z = 13.1826/(1 + DOC)^0.61 # log(z) = 1.12 - 0.61log(DOC+1)
## P mg/m2 phosphorus
## R mg/m2 recycled phosphorus
## A mg/m2 chlorophyll


function beisner!(du,u,h,p,t)
    A, P, B = u
    α, v, z, g, s, r, q, a, f, b, m, H, L, tau = p

    R = (r*(A^q))/(k^q + A^q)
    W = 3.2*exp(-3.6926*B)
    C = W < 2.58 ? 0.5 : 0.35 # mm

    Blag3 = h(p, t-tau)[1]

    du[1] = dA = α*P*A - ((v/z) + g*C + s)*A
    du[2] = dP = L + R - α*P*A - s*P
    du[3] = dB = B*exp(a*Blag3^f/(b^f+Blag3^f))-m*B-H*B
end

u0 = [1.0;0.0;1.5]
h(p,t) = u0[3] # assume B0 for initial delay values
tau = 3
lags = [tau]
p = (α, v, z, g, s, r, q, a, f, b, m, H, L, tau)
tspan = (0.0,10.0)
prob = DDEProblem(beisner!,u0,h,tspan,p;constant_lags=lags)
sol = solve(prob, MethodOfSteps(Tsit5()))

plot(sol.t, [a[1] for a in sol.u])
plot!(sol.t, [a[2] for a in sol.u])
gui()


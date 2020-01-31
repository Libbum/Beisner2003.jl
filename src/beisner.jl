using DifferentialEquations
using Plots

# Parameters
const f = 2;
const m = 0.01; # yr⁻¹
const a = 0.37;
const r = 8.7; # mg⋅m⁻²⋅d⁻¹
const q = 8; # Or 2
const k = 90; # mg⋅m⁻²
const s = 0.001; # d⁻¹
const b = 0.5;
const α = 0.8; # mg⋅m⁻²⋅d⁻¹
const v = 0.07; # m⋅d⁻¹
const g = 0.1; # d⁻¹


z = 0.04;
L = 1.429; # Exp(0.357)
DOC = 4.74; # Exp(1.556)
C = 0.35; # mc

function beisner!(du,u,p,t)
    A, P = u
    α, v, z, g, s, r, q, L, C = p
    R = (r*(A^q))/(k^q + A^q)
    du[1] = dA = α*P*A - ((v/z) + g*C + s)*A
    du[2] = dP = L + R - α*P*A - s*P
end

u0 = [1.0;0.0]
p = [α, v, z, g, s, r, q, L, C]
tspan = (0.0,10.0)
prob = ODEProblem(beisner!,u0,tspan,p)
sol = solve(prob)

plot(sol)
gui()

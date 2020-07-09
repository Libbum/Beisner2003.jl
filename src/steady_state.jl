using NLsolve
using Plots

# Parameters
const f = 2.0 # depensation in piscivores (B)
const m = 0.01 # yr⁻¹ natural piscivore mortality
const a = 0.37 # B recruitment assymptote (shape parameter)
const r = 8.7 # mg⋅m⁻²⋅d⁻¹ maximum recycling of P
const q = 8 # Or 2, steepness coefficient of R
const k = 90 # mg⋅m⁻² half saturation constant of P recycling #h in matlab
const s = 0.001 # d⁻¹ flushing rate # matlab q1
const b = 0.5 # adult density when recruitment = a/2 (B recruitment shape parameter)
const α = 0.8 # mg⋅m⁻²⋅d⁻¹ maximum growth rate of A per P input rate
const v = 0.07 # m⋅d⁻¹ sinking rate per algal cell
const g = 0.1 # d⁻¹ C grazing coefficient

# Drivers
# H yr⁻¹ harvest rate
#H = [0.1, 0.4]
# DOC mg/m2 dissolved organic carbon
DOC = [6.0, 6.0, 18.0, 18.0]
# L mg⋅m⁻²⋅d⁻¹ P loading
Ploading = 0.0:0.1:6.0

# State variables
# B g/m2 piscivore population biomass (pike)
## W g/m2 planktivore population biomass (bream)
#W = @. 3.2*exp(-3.6926*B)
## C mm zooplankton community length
#C = W < 2.58 ? 0.5 : 0.35 # mm
C = [0.5, 0.35, 0.5, 0.35] # For steady state, we assume a low or high harvest rate.
## Z m depth of thermocline
z = @. 13.1826/(1 + DOC)^0.61 # log10(z) = 1.12 - 0.61log10(DOC+1)
## P mg/m2 phosphorus
## R mg/m2 recycled phosphorus
## A mg/m2 chlorophyll

function beisnernl!(F, x, L, z, C)
    R = (r*(x[1]^q))/(k^q + x[1]^q)
    F[1] = α*x[2]*x[1] - ((v/z) + g*C + s)*x[1]
    F[2] = L + R - α*x[2]*x[1] - s*x[2]
end

handle = Vector{Plots.Plot}(undef,4)

for i in 1:4
    stable = Vector{Tuple{Float64,Float64}}(undef,0)
    unstable = Vector{Tuple{Float64,Float64}}(undef,0)
    for L in Ploading
        # Lower bound
        sol = nlsolve((F,x) -> beisnernl!(F,x, L, z[i], C[i]), [1.0;0.0])
        if converged(sol)
            push!(stable, (L, sol.zero[1]))
        end
        # Upper bound
        sol = nlsolve((F,x) -> beisnernl!(F,x, L, z[i], C[i]), [300.0;0.0])
        if converged(sol)
            push!(stable, (L, sol.zero[1]))
        end
        # Instability region
        sol = nlsolve((F,x) -> beisnernl!(F,x, L, z[i], C[i]), [70.0;0.0])
        if converged(sol)
            push!(unstable, (L, sol.zero[1]))
        end
    end
    handle[i] = scatter(stable; legend = false, markersize = 3)
    xlims!((0,6))
    ylims!((0, 300))
    scatter!(handle[i], unstable; legend = false, color = :red, markershape=:hline)
end

plot(handle..., layout = (2,2))


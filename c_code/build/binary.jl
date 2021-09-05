using DifferentialEquations
using Plots

#This cell contains the derivative of the Hamiltonian, used to calculate change in position and change in momentum

function Hamiltonian(du, u, p, t)
    x1 = u[1]
    y1 = u[2]
    p1x = u[3]
    p1y = u[4]
    x2 = u[5]
    y2 = u[6]
    p2x = u[7]
    p2y = u[8]

    m1 = p[1]
    m2 = p[2]
    G = p[3]

	difX, difY = x2 - x1, y2 - y1
	r = (((difX ^ 2) + (difY ^ 2)) ^ 1.5)
	top = G * m1 * m2

	du[1], du[2] = p1x / m1, p1y / m1
	du[3], du[4] = top * difX / r, top * difY / r

	du[5], du[6] = p2x / m2, p2y / m2
	du[7], du[8] = - top * difX / r, - top * difY / r
end

#Ensure command line arguments are present
if (size(ARGS, 1) != 7)
	println("Please provide first the masses of the initial bodies in solar masses, then the initial separation between them in kilometers, then eccentricity, then the y boundary for the plot, then the negative x boundary and the postivie x boundary.")
	exit()
end

#Read in the data from the input file
m1, m2, D, ecc, ymax, xmin, xmax = parse(Float64, ARGS[1]), parse(Float64, ARGS[2]), parse(Float64, ARGS[3]), parse(Float64, ARGS[4]), parse(Float64, ARGS[5]), parse(Float64, ARGS[6]), parse(Float64, ARGS[7])

G = 1; # Gravitational constant

GSI = 6.647e-11
MSUN = 1.989e30
KM = 1000

M = MSUN
L = KM
T = sqrt(L^3 / (GSI * M))

pUnits = M * L / T
KEUnits = pUnits * L / T

## Calculates factors based on eccentricity

if !(0 <= ecc <= 1)
    ecc = 0.0
end
qfac = 1 - ecc
pfac = sqrt((1 + ecc)/(1 - ecc))

mu = m1 * m2 / (m1 + m2)
p_theta = sqrt(G * m1 * m2 * mu * D )

x2 = D / (1 + (m1 / m2))
x1 = x2 - D
py2 = p_theta / D
py1 = -py2

c0 = [m1, m2, G]
q01 = [x1 * qfac, 0.0]
p01 = [0.0, py1 * pfac]
q02 = [x2 * qfac, 0.0]
p02 = [0.0, py2 * pfac]
u0 = collect(Base.Iterators.flatten([q01, p01, q02, p02]))

pi1, pi2 = sqrt((p01[1]^2) + (p01[2]^2)), sqrt((p02[1]^2) + (p02[2]^2))
KEi1, KEi2 = pi1^2 / (2 * m1), pi2^2 / (2 * m2)
initialKE1, initialKE2 = KEi1 * KEUnits, KEi2 * KEUnits
println("Initial kinetic energies, in joules, are ", initialKE1, " and ", initialKE2)

tspan = (0.0, 500.0)
prob = ODEProblem(Hamiltonian, u0, tspan, c0)
sol = solve(prob, RK4(), adaptive = true, reltol = 1.0e-8, abstol = 1.0e-8, saveat = 0.1)
ss1 = sol[1,:]
ss2 = sol[2,:]
ss3 = sol[5,:]
ss4 = sol[6,:]

len = size(sol[1,:])[1]

pf1, pf2 = sqrt((sol[3,:][len]^2) + (sol[4,:][len]^2)), sqrt((sol[7,:][len]^2) + (sol[8,:][len]^2))
KEf1, KEf2 = pf1^2 / (2 * m1), pf2^2 / (2 * m2)
finalKE1, finalKE2 = KEf1 * KEUnits, KEf2 * KEUnits
dKE1, dKE2 = finalKE1 - initialKE1, finalKE2 - initialKE2
xf, yf = ss1[len] - ss3[len], ss2[len] - ss4[len]
Df = sqrt(xf^2 + yf^2)
println("Final kinetic energies, in joules, are ", finalKE1, " and ", finalKE1)
println("Inital and final distances between bodies, in kilometers, are ", D, " and ", Df)

@userplot TwoBodyPlot
@recipe function f(tb::TwoBodyPlot)
    x, y, xarray, yarray, negx, posx, negy, posy = tb.args
    n = 100
    xlims --> (negx, posx)
    ylims --> (negy, posy)
    append!(xarray, x)
    append!(yarray, y)
    if size(xarray)[1] > n
        deleteat!(xarray, 1)
        deleteat!(yarray, 1)
    end
    linewidth --> range(0, 10, length = n)
    #seriesalpha --> range(0, 1, length = n)
    aspect_ratio --> 1
    label --> false
    xarray, yarray
end

xarr1, yarr1, xarr2, yarr2 = [], [], [], []

anim = @gif for i = 1:size(ss1)[1]
    twobodyplot(ss1[i], ss2[i], xarr1, yarr1, xmin, xmax, -ymax, ymax)
    twobodyplot!(ss3[i], ss4[i], xarr2, yarr2, xmin, xmax, -ymax, ymax)
end every 5

display(anim)

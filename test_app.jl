using DifferentialEquations
using DelimitedFiles
using Pkg
Pkg.activate("FormattedEquations")
using FormattedEquations

Base.@ccallable function julia_main()::Cint
	try
		real_main()
	catch
		Base.invokelatest(Base.display_error, Base.catch_stack())
		return 1
	end
	return 0
end

function real_main()::Cint
	#=
	file = ARGS[1]
	arr = readdlm(file, ' ', Float64, '\n')
	G, M, L, T = arr[1,1], arr[1,2], arr[1,3], arr[1,4]
	=#

	file = ARGS[1]
	arr = readdlm(file, ' ', Float64, '\n')
	nbody = size(arr)[1] - 1
	G, M, L, T = arr[1,1], arr[1,2], arr[1,3], arr[1,4]
	data = arr[2,:][2:end]
	for i in 2:nbody
		append!(data, arr[i+1,:][2:end])
	end
	c0 = arr[:,1][2:end]
	append!(c0,G)

	m1, x1, y1, z1, px1, py1, pz1 = arr[2,:]
	m2, x2, y2, z2, px2, py2, pz2 = arr[3,:]

    C = 1.000000000000; #Speed of light
    tspan = (0.0, 1 * 12000.0); # The amount of time for which the simulation runs
	xi, yi, zi = x1 - x2, y1 - y2, z1 - z2
	D = sqrt(xi^2 + yi^2 + zi^2)

	h01 = [0]

	#TAYDEN WUZ HERE

	#=
	schwarz = [R_Schwarz(C,G,m1), R_Schwarz(C,G,m2)]

	condition_collision(u,t,integrator) = schwarz[1] + schwarz[2] >= sqrt((u[1] - u[5])^2 + (u[2] - u[6])^2)
	affect_collision!(integrator) = terminate!(integrator)

	cb_collision = DiscreteCallback(condition_collision, affect_collision!, save_positions=(true,true))
	=# #FIXME Make this section with collision detection n-body compatible.

	#u0 = collect(Base.Iterators.flatten([q01, p01, q02, p02, h01]));
	u0 = collect(Base.Iterators.flatten([data, h01]));

    prob = ODEProblem(FormattedEquations.PM, u0, tspan, c0);
    sol = DifferentialEquations.solve(prob, Vern9(), #=callback=cb_collision,=# reltol = 1.0e-9, abstol = 1.0e-9, saveat = 10);

    #=

    xlist1 = sol[1,:]
    ylist1 = sol[2,:]
    zlist1 = sol[3,:]
    xlist2 = sol[7,:]
    ylist2 = sol[8,:]
    zlist2 = sol[9,:]

    momentumx1 = sol[4,:]
    momentumy1 = sol[5,:]
    momentumz1 = sol[6,:]
    momentumx2 = sol[10,:]
    momentumy2 = sol[11,:]
    momentumz2 = sol[12,:]

    hamilArr = sol[(nbody * 6) + 1,:];

    deleteat!(hamilArr, 1);
    originalHamiltonian = copy(hamilArr[1])
    lenHamil = size(hamilArr)
    originalArr = ones(lenHamil) * originalHamiltonian
    hamilVariance = hamilArr - originalArr

    xdis, ydis, zdis = xlist1 - xlist2, ylist1 - ylist2, zlist1 - zlist2
    distArr = (xdis.^2 + ydis.^2 + zdis.^2).^0.5
    distVar = distArr - ( ones(size(distArr)) * D )

    linear_momentum1 = ((copy(momentumx1) .*= momentumx1) + (copy(momentumy1) .*= momentumy1) + (copy(momentumz1) .*= momentumz1)).^0.5
    linear_momentum2 = ((copy(momentumx2) .*= momentumx2) + (copy(momentumy2) .*= momentumy2) + (copy(momentumz2) .*= momentumz2)).^0.5

    angular_momentum1 = copy(distArr) .*= linear_momentum1
    angular_momentum2 = copy(distArr) .*= linear_momentum2

    xf, yf, zf = xlist1[end] - xlist2[end], ylist1[end] - ylist2[end], zlist1[end] - zlist2[end]
    Df = sqrt(xf^2 + yf^2 + zf^2)
    println("Final distance between objects = ", Df)
	#println("Schwarzchild distance is ", schwarz[1] + schwarz[2])
    println()

	=#

	open("PMdata.csv", "w+") do io
		writedlm(io, sol, ',')
	end

    return 0
end

condition_orbits(u,t,integrator) = (u[1] == 0) && (u[2] > 0)
affect_orbits!(integrator) = integrator.u[10] += 1

CSI = 3.00e8;
GSI = 6.647e-11;
MSUN = 1.989e30;

function R_Schwarz(C, G, M)
	return 2 * G * M / (C^2)
end

julia_main()

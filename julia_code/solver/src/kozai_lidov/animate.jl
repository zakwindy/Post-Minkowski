using Plots
using CSV
using DataFrames

@userplot TwoBodyPlot
@recipe function f(tb::TwoBodyPlot)
    m, x, y, z, xarray, yarray, zarray, negx, posx, negy, posy, negz, posz, r = tb.args
    n = r
    xlims --> (negx, posx)
    ylims --> (negy, posy)
    zlims --> (negz, posz)
    append!(xarray, x)
    append!(yarray, y)
    append!(zarray, z)
    if size(xarray)[1] > n
        deleteat!(xarray, 1)
        deleteat!(yarray, 1)
        deleteat!(zarray, 1)
    end
    linewidth --> range(0, m * 5, length = n)
    #seriesalpha --> range(0, 1, length = n)
    aspect_ratio --> 1
    label --> false
    xarray, yarray, zarray
end

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
	m1 = 0.0;
	m2 = 0.0;
	m3 = 0.0;
	try 			# make sure that a file has been included
		f = open(ARGS[1])
		m1 = ARGS[2]
		m2 = ARGS[3]
		m3 = ARGS[4]
		close(f)
	catch e
		println("Please include a .csv file output by the solver as the first argument.")
		return 1
	end

	marr = [m1, m2, m3]
	m0arr = [0.5, 0.5, 0.5]

	if !((m1 == m2) && (m2 == m3))
		m0arr[findmax(marr)[2]] = 1
		m0arr[findmin(marr)[2]] = 0.25
	end

	file = ARGS[1]
	arr = CSV.read(file, DataFrame)

	C = 1.0			# set relativistic C
	mSun = 1.0		# set the mass to be in solar mass units

	C_CGS = 2.998e10;
	G_CGS = 6.674e-8;
	mSun_CGS = 1.989e33;
	AU_CGS = 1.496e14; # 1 AU in cm

	G = 1.0;
	M = mSun_CGS;		#units of mass
	L = M * (G_CGS / G) * ((C / C_CGS)^2);		#units of length
	T = L * C / C_CGS;		#units of time

	qx1 = arr[:,"qx1"]*L/AU_CGS
	qy1 = arr[:,"qy1"]*L/AU_CGS
	qz1 = arr[:,"qz1"]*L/AU_CGS
	qx2 = arr[:,"qx2"]*L/AU_CGS
	qy2 = arr[:,"qy2"]*L/AU_CGS
	qz2 = arr[:,"qz2"]*L/AU_CGS
	qx3 = arr[:,"qx3"]*L/AU_CGS
	qy3 = arr[:,"qy3"]*L/AU_CGS
	qz3 = arr[:,"qz3"]*L/AU_CGS

	maxx = [findmax(qx1)[1],findmax(qx2)[1],findmax(qx3)[1]]
	minx = [findmin(qx1)[1],findmin(qx2)[1],findmin(qx3)[1]]
	maxy = [findmax(qy1)[1],findmax(qy2)[1],findmax(qy3)[1]]
	miny = [findmin(qy1)[1],findmin(qy2)[1],findmin(qy3)[1]]
	maxz = [findmax(qz1)[1],findmax(qz2)[1],findmax(qz3)[1]]
	minz = [findmin(qz1)[1],findmin(qz2)[1],findmin(qz3)[1]]

	maxX = findmax(maxx)[1]
	minX = findmin(minx)[1]
	maxY = findmax(maxy)[1]
	minY = findmin(miny)[1]
	maxZ = findmax(maxz)[1]
	minZ = findmin(minz)[1]

	xarr1, yarr1, zarr1, xarr2, yarr2, zarr2, xarr3, yarr3, zarr3 = [], [], [], [], [], [], [], [], []

	anim = @animate for i = 1:size(qx1)[1]
    	twobodyplot(m0arr[1], qx1[i], qy1[i], qz1[i], xarr1, yarr1, zarr1, minX, maxX, minY, maxY, minZ, maxZ, 100)
    	twobodyplot!(m0arr[2], qx2[i], qy2[i], qz2[i], xarr2, yarr2, zarr2, minX, maxX, minY, maxY, minZ, maxZ, 100)
		twobodyplot!(m0arr[3], qx3[i], qy3[i], qz3[i], xarr3, yarr3, zarr3, minX, maxX, minY, maxY, minZ, maxZ, 100)
	end every 10

	gif(anim, "animation.gif", fps = 15)

	return 0
end

julia_main()

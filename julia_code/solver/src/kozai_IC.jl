using DelimitedFiles

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
	G = 1.0			# set relativistic G
	C = 1.0			# set relativistic C
	mSun = 1.0		# set the mass to be in solar mass units

	try 			# make sure that a file has been included
		f = open(ARGS[1])
		close(f)
	catch e
		println("Please include a file with the masses of the three bodies as the first line (the inner orbit bodies first and the outer body third, followed by three zeros to keep the dimensions of the array constant) and the orbital parameters as the second line and third lines, first semi-major axis, then eccentricity, then inclination, then argument of periastron, then longitude of ascending node, and then the mean anomaly.")
		return 1
	end

	file = ARGS[1]
	arr = readdlm(file, ' ', Float64, '\n')

	rows = size(arr)[1]			#make sure the dimensions of the file are correct
	columns = size(arr)[2]
	if rows != 3
		println("Incorrect dimensions in initial file, needs 3 rows.")
		return 2
	end
	if columns != 6
		println("Incorrect dimensions in initial file, needs 6 columns.")
		return 2
	end

	m1, m2, m3 = arr[1,1], arr[1,2], arr[1,3]		# set the masses
	aIn, eIn, iIn, wIn, omegaIn, mIn = arr[2,1], arr[2,2], arr[2,3], arr[2,4], arr[2,5], arr[2,6]		# set the inner orbit parameters
	aOut, eOut, iOut, wOut, omegaOut, mOut = arr[3,1], arr[3,2], arr[3,3], arr[3,4], arr[3,5], arr[3,6]		# set the outer orbit parameters

	return 0
end

function initialdata(m1,m2,userG,a,e,i,w,omega,M,tol0)

    mu=m1*m2/(m1+m2)

    u0 = 10;        # initializing u0
    tol  = 1;      # Initalizing iterative tolerance

    while abs(tol) > tol0
        u0 = u0 - (-M + u0 - e*sin(u0))/(1 - e*cos(u0))
        tol = -M + u0 - e*sin(u0)
    end

    # Using eccentric anomaly to find true anomaly

    f = atan(sin(u0)*sqrt(1-e^2)/(cos(u0)-e));

    # Solving for Polar Coordinates

    r     = a*(1-e^2)/(1-e*cos(f));
    phi   = omega + atan(tan(w+f)*cos(i));
    theta = acos(sin(w+f)*sin(i));

    # Some useful quantities in determining velocities

    gr     = a*(1-e^2)*e*sin(f)/(1+e*cos(f))^2;
    gtheta = (-1/sin(theta))*cos(w+f)*sin(i);
    gphi   = cos(phi-omega)^2*cos(i)/cos(w+f)^2;

    fdot   = sqrt(userG*mbi*((2/r)-(1/a))*1/(gr^2 + (r*gtheta)^2 + (r*sin(theta)*gphi)^2))

    # Setting up polar velocities

    rdot     = gr*fdot;
    thetadot = gtheta*fdot;
    phidot   = gphi*fdot;

    ### TODO: Convert Polar speeds and coordinates to cartesian
    #         Determine where exactly this function can be used
    #         Implement this function to get initial data
    #         Return correct data points
    #         Add conserved quantities
    ###

    x=r*sin(theta)*cos(phi)
    y=r*sin(theta)*sin(phi)
    z=r*cos(theta)

    # Double check math for everything past this point

    xdot = rdot*sin(theta)*cos(phi) + r*thetadot*cos(theta)*cos(phi) - phidot*r*sin(phi)*sin(theta)
    ydot = rdot*sin(theta)*sin(phi) + r*thetadot*cos(theta)*sin(phi) + r*phidot*sin(theta)*cos(phi)
    zdot = rdot*cos(theta) - r*thetadot*sin(theta)

    x1 = m2*x/(m1+m2)
    y1 = m2*y/(m1+m2)
    z1 = m2*z/(m1+m2)

    x2 = -m1*x/(m1+m2)
    y2 = -m1*y/(m1+m2)
    z2 = -m1*z/(m1+m2)

    px1 = mu*xdot
    py1 = mu*ydot
    pz1 = mu*zdot

    px2 = -mu*xdot
    py2 = -mu*ydot
    pz2 = -mu*zdot

    [[x1,y1,z1],[x2,y2,z2],[px1,py1,pz1],[px2,py2,pz2]]



    # Calculate angular momentum
    # linear momentum = mu * v
    # mv1 + mv2 = 0
    # Consider that c=1, G=1
    # Double check that using SI unites

end

julia_main()

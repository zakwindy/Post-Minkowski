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
	C = 1.0			# set relativistic C
	mSun = 1.0		# set the mass to be in solar mass units

	C_CGS = 2.998e10;
	G_CGS = 6.674e-8;
	mSun_CGS = 1.989e33;
	AU = 1.496e14; # 1 AU in cm

	## Here the code units are defined in terms of CGS units


	# Defined using C = G = 1
	G = 1.0;
	M = mSun_CGS;		#units of mass
	L = M * (G_CGS / G) * ((C / C_CGS)^2);		#units of length
	T = L * C / C_CGS;		#units of time

	#=
	# Defined using T = 1 year
	M = mSun_CGS;
	T = 3600*24*365;
	L = T * C_CGS / C;
	G = M*G_CGS*((C/C_CGS)^2)/L;
	=#
	#=
	# Defined using L
	M = mSun_CGS;
	L = 1 * AU;
	T = L * C / C_CGS;
	G = M*G_CGS*((C/C_CGS)^2)/L;
	=#



	println("M = ", M);
	println("L = ", L);
	println("T = ", T);
	println("G = ", G);

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
	aInAU, eIn, iIn_degree, wIn_degree, omegaIn, MIn_degree = arr[2,1], arr[2,2], arr[2,3], arr[2,4], arr[2,5], arr[2,6]		# set the inner orbit parameters
	aOutAU, eOut, iOut_degree, wOut_degree, omegaOut, MOut_degree = arr[3,1], arr[3,2], arr[3,3], arr[3,4], arr[3,5], arr[3,6]		# set the outer orbit parameters
	tol0 = 1e-30;		# set tolerance level
	mbi = m1 + m2;

	aIn = aInAU * AU / L;		#convert from AU to code units
	aOut = aOutAU * AU / L;
	iIn = iIn_degree * pi / 180;	#convert from degrees to radians
	iOut = iOut_degree * pi / 180;
	wIn = wIn_degree * pi / 180;
	wOut = wOut_degree * pi / 180;
	MIn = MIn_degree * pi / 180;
	MOut = MOut_degree * pi / 180;

	valuesIn=initialdata(m1,m2,mbi,G,aIn,eIn,iIn,wIn,omegaIn,MIn,tol0);
	valuesOut=initialdata(m3,mbi,mbi,G,aOut,eOut,iOut,wOut,omegaOut,MOut,tol0);

	r1 = valuesIn[1]+valuesOut[2]
	r2 = valuesIn[2]+valuesOut[2]
	r3 = valuesOut[1]
	rdot1 = valuesIn[3]+valuesOut[4]
	rdot2 = valuesIn[4]+valuesOut[4]
	rdot3 = valuesOut[3]

	p1 = m1*rdot1
	p2 = m2*rdot2
	p3 = m3*rdot3

	output = open("ICfile0", "w+")

	write(output, string(G, ' ', M, ' ', L, ' ', T,  " 0 0 0\n"));
	write(output, string(m1,' ',r1[1],' ',r1[2],' ',r1[3],' ',p1[1],' ',p1[2],' ',p1[3],'\n'))
	write(output, string(m2,' ',r2[1],' ',r2[2],' ',r2[3],' ',p2[1],' ',p2[2],' ',p2[3],'\n'))
	write(output, string(m3,' ',r3[1],' ',r3[2],' ',r3[3],' ',p3[1],' ',p3[2],' ',p3[3],'\n'))

	return 0
end

function initialdata(m1,m2,mBin,userG,a,e,i,w,omega,M,tol0)

    mu=m1*m2/(m1+m2)
		Mtot = m1+m2

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

    gr     = a*(1-e^2)*e*sin(f)/((1+e*cos(f))^2);
    gtheta = (-1/sin(theta))*cos(w+f)*sin(i);
    gphi   = cos(phi-omega)^2*cos(i)/(cos(w+f)^2);

    fdot   = sqrt(userG*mBin*((2/r)-(1/a))*1/(gr^2 + (r*gtheta)^2 + (r*sin(theta)*gphi)^2))

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

	Lx = mu*(y*zdot - z*ydot);
	Ly = mu*(x*zdot - z*xdot);
	Lz = mu*(x*ydot - y*xdot);

    x1 = m2*x/(m1+m2)
    y1 = m2*y/(m1+m2)
    z1 = m2*z/(m1+m2)

    x2 = -m1*x/(m1+m2)
    y2 = -m1*y/(m1+m2)
    z2 = -m1*z/(m1+m2)

    px1 = 1*mu*xdot
    py1 = 1*mu*ydot
    pz1 = 1*mu*zdot

    px2 = -px1
    py2 = -py1
    pz2 = -pz1

	x1dot = px1/m1
	y1dot = py1/m1
	z1dot = pz1/m1

	x2dot = px2/m2
	y2dot = py2/m2
	z2dot = pz2/m2

    [[x1,y1,z1],[x2,y2,z2],[x1dot,y1dot,z1dot],[x2dot,y2dot,z2dot]]



    # Calculate angular momentum
    # linear momentum = mu * v
    # mv1 + mv2 = 0
    # Consider that c=1, G=1
    # Double check that using SI unites

end

julia_main()

using DifferentialEquations
using LSODA
using DelimitedFiles

function julia_main()::Cint
	#=
	file = ARGS[1]
	arr = readdlm(file, ' ', Float64, '\n')
	G, M, L, T = arr[1,1], arr[1,2], arr[1,3], arr[1,4]
	=#

	file = ARGS[1]
	arr = readdlm(file, ' ', Float64, '\n')
	nbody = size(arr)[1]
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

    #c0 = [m1, m2, G]

	q01 = [x1, y1]
    p01 = [px1, py1]
    q02 = [x2, y2]
    p02 = [px2, py2]

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

    prob = ODEProblem(PM, u0, tspan, c0);
    sol = DifferentialEquations.solve(prob, Vern9(), #=callback=cb_collision,=# reltol = 1.0e-9, abstol = 1.0e-9, saveat = 10);

    len = size(sol[1,:])[1]

    xlist1 = sol[1,:]
    ylist1 = sol[2,:]
    zlist1 = zeros(len)
    xlist2 = sol[5,:]
    ylist2 = sol[6,:]
    zlist2 = zeros(len)

    momentumx1 = sol[3,:]
    momentumy1 = sol[4,:]
    momentumz1 = zeros(len)
    momentumx2 = sol[7,:]
    momentumy2 = sol[8,:]
    momentumz2 = zeros(len)

    hamilArr = sol[9,:];

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

    xf, yf, zf = xlist1[len] - xlist2[len], ylist1[len] - ylist2[len], zlist1[len] - zlist2[len]
    Df = sqrt(xf^2 + yf^2 + zf^2)
    println("Final distance between objects = ", Df)
	println("Schwarzchild distance is ", schwarz[1] + schwarz[2])
    println()

	open("PMdata.csv", "w+") do io
		writedlm(io, [xlist1, ylist1, zlist1, xlist2, ylist2, zlist2, hamilVariance, distVar, linear_momentum1, linear_momentum2, angular_momentum1, angular_momentum2], ',')
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

function PM(du, u, p, t)
	qax = u[1]
	qay = u[2]
	pax = u[3]
	pay = u[4]
	qbx = u[5]
	qby = u[6]
	pbx = u[7]
	pby = u[8]

	ma = p[1]
	mb = p[2]
	G = p[3]

	#Optimization variables
	o3=ma*ma
	o4=pax*pax
	o5=pay*pay
	o6=o3+o4+o5
	o7=1/sqrt(o6)
	o16=mb*mb
	o17=pbx*pbx
	o18=pby*pby
	o19=o16+o17+o18
	o20=sqrt(o19)
	o10=o4+o5
	o13=1/o6
	o21=-qbx
	o22=o21+qax
	o23=o22*o22
	o24=-qby
	o25=o24+qay
	o26=o25*o25
	o27=o23+o26
	o28=1/sqrt(o27)
	o9=sqrt(o6)
	o11=o6^-2
	o12=-2*o10*o11*pax
	o14=2*o13*pax
	o15=o12+o14
	o30=o10*o13
	o31=o17+o18
	o32=1/o19
	o33=o31*o32
	o34=1+o30+o33
	o36=-qax
	o37=o36+qbx
	o38=o37*o37
	o39=-qay
	o40=o39+qby
	o41=o40*o40
	o42=o38+o41
	o43=1/sqrt(o42)
	o48=7*pbx
	o73=pax*pbx
	o74=pay*pby
	o75=o73+o74
	o76=o75*o75
	o64=o22*o28*pax
	o65=o25*o28*pay
	o66=o64+o65
	o67=o66*o66
	o49=o22*o28*pbx
	o50=o25*o28*pby
	o51=o49+o50
	o68=o3+o67
	o69=sqrt(o68)
	o81=o51*o51
	o89=o67*o81
	o101=1/sqrt(o68)
	o63=1/sqrt(o19)
	o70=o69*o7
	o71=1+o70
	o77=-(o10*o76)
	o78=2*o67*o76
	o79=-2*o10*o51*o66*o75
	o80=o10*o10
	o82=o80*o81
	o83=o77+o78+o79+o82
	o102=o6^-1.5
	o85=o10*o31
	o86=-3*o31*o67
	o87=8*o51*o66*o75
	o88=-3*o10*o81
	o90=o85+o86+o87+o88+o89
	o92=-o17
	o93=-o18
	o94=o92+o93
	o125=2*o22*o28*o66*o81
	o107=o71^-2
	o84=2*o13*o83
	o91=o69*o7*o90
	o95=o67*o94
	o96=2*o51*o66*o75
	o97=-(o10*o81)
	o98=o76+o89+o95+o96+o97
	o99=2*o98
	o100=o84+o91+o99
	o55=o37*o43*pbx
	o56=o40*o43*pby
	o57=o55+o56
	o140=o57*o57
	o141=o140+o16
	o149=o37*o43*pax
	o150=o40*o43*pay
	o151=o149+o150
	o143=sqrt(o141)
	o128=2*o31*pax
	o120=2*o75*pbx
	o162=2*o140*o151*o37*o43
	o142=1/sqrt(o141)
	o144=o143*o63
	o145=1+o144
	o146=o145^-2
	o148=o31*o31
	o174=o151*o151
	o186=o140*o174
	o200=-2*o10*o11*pay
	o201=2*o13*pay
	o202=o200+o201
	o209=7*pby
	o72=o71^-3
	o239=2*o25*o28*o66*o81
	o138=o68^-1.5
	o242=2*o31*pay
	o234=2*o75*pby
	o264=2*o140*o151*o40*o43
	o173=-(o31*o76)
	o175=o148*o174
	o176=-2*o151*o31*o57*o75
	o177=2*o140*o76
	o178=o173+o175+o176+o177
	o179=2*o178*o32
	o180=-(o174*o31)
	o181=2*o151*o57*o75
	o182=-o4
	o183=-o5
	o184=o182+o183
	o185=o140*o184
	o187=o180+o181+o185+o186+o76
	o188=2*o187
	o189=-3*o174*o31
	o190=8*o151*o57*o75
	o191=-3*o10*o140
	o192=o186+o189+o190+o191+o85
	o193=o143*o192*o63
	o194=o179+o188+o193
	o280=o19^-2
	o281=-2*o280*o31*pbx
	o282=2*o32*pbx
	o283=o281+o282
	o290=7*pax
	o313=2*o22*o28*o51*o67
	o316=2*o75*pax
	o308=2*o10*pbx
	o329=2*o174*o37*o43*o57
	o299=o19^-1.5
	o364=-2*o280*o31*pby
	o365=2*o32*pby
	o366=o364+o365
	o373=7*pay
	o395=2*o25*o28*o51*o67
	o398=2*o75*pay
	o390=2*o10*pby
	o411=2*o174*o40*o43*o57
	o356=o145^-3
	o358=o141^-1.5
	o443=o27^-1.5
	o445=o42^-1.5
	o449=7*o75
	o458=-(o23*o443*pax)
	o459=o28*pax
	o460=-(o22*o25*o443*pay)
	o461=o458+o459+o460
	o453=-(o23*o443*pbx)
	o454=o28*pbx
	o455=-(o22*o25*o443*pby)
	o456=o453+o454+o455
	o493=2*o456*o51*o67
	o494=2*o461*o66*o81
	o470=o38*o445*pax
	o471=o37*o40*o445*pay
	o472=-(o43*pax)
	o473=o470+o471+o472
	o465=o38*o445*pbx
	o466=o37*o40*o445*pby
	o467=-(o43*pbx)
	o468=o465+o466+o467
	o519=2*o174*o468*o57
	o520=2*o140*o151*o473
	o450=o51*o66
	o451=o449+o450
	o477=o151*o57
	o478=o449+o477
	o543=o28*pay
	o544=-(o22*o25*o443*pax)
	o545=-(o26*o443*pay)
	o546=o543+o544+o545
	o548=o28*pby
	o549=-(o22*o25*o443*pbx)
	o550=-(o26*o443*pby)
	o551=o548+o549+o550
	o580=2*o546*o66*o81
	o583=2*o51*o551*o67
	o507=1/o68
	o561=o37*o40*o445*pax
	o562=o41*o445*pay
	o563=-(o43*pay)
	o564=o561+o562+o563
	o556=o37*o40*o445*pbx
	o557=o41*o445*pby
	o558=-(o43*pby)
	o559=o556+o557+o558
	o607=2*o174*o559*o57
	o608=2*o140*o151*o564
	o532=1/o141
	o636=o23*o443*pax
	o637=-(o28*pax)
	o638=o22*o25*o443*pay
	o639=o636+o637+o638
	o631=o23*o443*pbx
	o632=-(o28*pbx)
	o633=o22*o25*o443*pby
	o634=o631+o632+o633
	o669=2*o51*o634*o67
	o670=2*o639*o66*o81
	o648=-(o38*o445*pax)
	o649=-(o37*o40*o445*pay)
	o650=o43*pax
	o651=o648+o649+o650
	o643=-(o38*o445*pbx)
	o644=-(o37*o40*o445*pby)
	o645=o43*pbx
	o646=o643+o644+o645
	o694=2*o174*o57*o646
	o695=2*o140*o151*o651
	o717=-(o28*pay)
	o718=o22*o25*o443*pax
	o719=o26*o443*pay
	o720=o717+o718+o719
	o722=-(o28*pby)
	o723=o22*o25*o443*pbx
	o724=o26*o443*pby
	o725=o722+o723+o724
	o754=2*o66*o720*o81
	o757=2*o51*o67*o725
	o735=-(o37*o40*o445*pax)
	o736=-(o41*o445*pay)
	o737=o43*pay
	o738=o735+o736+o737
	o730=-(o37*o40*o445*pbx)
	o731=-(o41*o445*pby)
	o732=o43*pby
	o733=o730+o731+o732
	o781=2*o174*o57*o733
	o782=2*o140*o151*o738

	#Equations for changes in body 1
	du[1] = 0.25*G*(o28*(o48+o22*o28*o51)+o43*(o48+o37*o43*o57))+o7*pax-0.5*G*(o15*o20*o28*o9+o15*o20*o43*o9+o20*o28*o34*o7*pax+o20*o34*o43*o7*pax)-0.25*G*(-(o100*o107*o138*o22*o28*o43*o63*o66)-o102*o142*o146*o194*o28*pax-2*o100*o101*o43*o63*o72*(o101*o22*o28*o66*o7-o102*o69*pax)+o142*o146*o28*o7*(2*(o120+o162-2*o151*o31*o37*o43+2*o37*o43*o57*o75-2*o140*pax+2*o151*o57*pbx)+o143*o63*(o128+o162-6*o151*o31*o37*o43+8*o37*o43*o57*o75-6*o140*pax+8*o151*o57*pbx)+2*o32*(2*o148*o151*o37*o43-2*o31*o37*o43*o57*o75-2*o151*o31*o57*pbx+4*o140*o75*pbx-2*o31*o75*pbx))+o101*o107*o43*o63*(o101*o22*o28*o66*o7*o90-4*o11*o83*pax-o102*o69*o90*pax+2*(o120+o125+2*o22*o28*o51*o75+2*o22*o28*o66*o94-2*o81*pax+2*o51*o66*pbx)+o69*o7*(o125+o128-6*o22*o28*o31*o66+8*o22*o28*o51*o75-6*o81*pax+8*o51*o66*pbx)+2*o13*(-2*o10*o22*o28*o51*o75+4*o22*o28*o66*o76-4*o51*o66*o75*pax-2*o76*pax+4*o10*o81*pax-2*o10*o51*o66*pbx-2*o10*o75*pbx+4*o67*o75*pbx)))

	du[2] = 0.25*G*(o28*(o209+o25*o28*o51)+o43*(o209+o40*o43*o57))+o7*pay-0.5*G*(o20*o202*o28*o9+o20*o202*o43*o9+o20*o28*o34*o7*pay+o20*o34*o43*o7*pay)-0.25*G*(-(o100*o107*o138*o25*o28*o43*o63*o66)-o102*o142*o146*o194*o28*pay-2*o100*o101*o43*o63*o72*(o101*o25*o28*o66*o7-o102*o69*pay)+o142*o146*o28*o7*(2*(o234+o264-2*o151*o31*o40*o43+2*o40*o43*o57*o75-2*o140*pay+2*o151*o57*pby)+o143*o63*(o242+o264-6*o151*o31*o40*o43+8*o40*o43*o57*o75-6*o140*pay+8*o151*o57*pby)+2*o32*(2*o148*o151*o40*o43-2*o31*o40*o43*o57*o75-2*o151*o31*o57*pby+4*o140*o75*pby-2*o31*o75*pby))+o101*o107*o43*o63*(o101*o25*o28*o66*o7*o90-4*o11*o83*pay-o102*o69*o90*pay+2*(o234+o239+2*o25*o28*o51*o75+2*o25*o28*o66*o94-2*o81*pay+2*o51*o66*pby)+o69*o7*(o239+o242-6*o25*o28*o31*o66+8*o25*o28*o51*o75-6*o81*pay+8*o51*o66*pby)+2*o13*(-2*o10*o25*o28*o51*o75+4*o25*o28*o66*o76-4*o51*o66*o75*pay-2*o76*pay+4*o10*o81*pay-2*o10*o51*o66*pby-2*o10*o75*pby+4*o67*o75*pby)))

	du[3] = -0.25*G*(-(o22*o443*o451)+o37*o445*o478+o43*(o151*o468+o473*o57)+o28*(o461*o51+o456*o66))+0.5*G*(-(o20*o22*o34*o443*o9)+o20*o34*o37*o445*o9)+0.25*G*(o100*o101*o107*o37*o445*o63-o100*o107*o138*o43*o461*o63*o66-o142*o146*o194*o22*o443*o7-o146*o194*o28*o358*o468*o57*o7-2*o194*o28*o356*o468*o532*o57*o63*o7-2*o100*o43*o461*o507*o63*o66*o7*o72+o142*o146*o28*o7*(o142*o192*o468*o57*o63+2*(-2*o151*o31*o473+o519+o520+2*o184*o468*o57+2*o151*o468*o75+2*o473*o57*o75)+o143*o63*(-6*o151*o31*o473+o519+o520-6*o10*o468*o57+8*o151*o468*o75+8*o473*o57*o75)+2*o32*(2*o148*o151*o473-2*o151*o31*o468*o75-2*o31*o473*o57*o75+4*o468*o57*o76))+o101*o107*o43*o63*(o69*o7*(o493+o494-6*o10*o456*o51-6*o31*o461*o66+8*o461*o51*o75+8*o456*o66*o75)+2*o13*(-2*o10*o461*o51*o75-2*o10*o456*o66*o75+4*o461*o66*o76+2*o456*o51*o80)+o101*o461*o66*o7*o90+2*(o493+o494-2*o10*o456*o51+2*o461*o51*o75+2*o456*o66*o75+2*o461*o66*o94)))

	du[4] = -0.25*G*(-(o25*o443*o451)+o40*o445*o478+o43*(o151*o559+o564*o57)+o28*(o51*o546+o551*o66))+0.5*G*(-(o20*o25*o34*o443*o9)+o20*o34*o40*o445*o9)+0.25*G*(o100*o101*o107*o40*o445*o63-o100*o107*o138*o43*o546*o63*o66-o142*o146*o194*o25*o443*o7-o146*o194*o28*o358*o559*o57*o7-2*o194*o28*o356*o532*o559*o57*o63*o7-2*o100*o43*o507*o546*o63*o66*o7*o72+o142*o146*o28*o7*(o142*o192*o559*o57*o63+2*(-2*o151*o31*o564+2*o184*o559*o57+o607+o608+2*o151*o559*o75+2*o564*o57*o75)+o143*o63*(-6*o151*o31*o564-6*o10*o559*o57+o607+o608+8*o151*o559*o75+8*o564*o57*o75)+2*o32*(2*o148*o151*o564-2*o151*o31*o559*o75-2*o31*o564*o57*o75+4*o559*o57*o76))+o101*o107*o43*o63*(o69*o7*(-6*o10*o51*o551+o580+o583-6*o31*o546*o66+8*o51*o546*o75+8*o551*o66*o75)+2*o13*(-2*o10*o51*o546*o75-2*o10*o551*o66*o75+4*o546*o66*o76+2*o51*o551*o80)+o101*o546*o66*o7*o90+2*(-2*o10*o51*o551+o580+o583+2*o51*o546*o75+2*o551*o66*o75+2*o546*o66*o94)))

	#Changes in body 2
	du[5] = 0.25*G*(o43*(o290+o151*o37*o43)+o28*(o290+o22*o28*o66))+o63*pbx-0.5*G*(o20*o28*o283*o9+o20*o283*o43*o9+o28*o34*o63*o9*pbx+o34*o43*o63*o9*pbx)-0.25*G*(-(o146*o194*o28*o358*o37*o43*o57*o7)-o100*o101*o107*o299*o43*pbx-2*o142*o194*o28*o356*o7*(o142*o37*o43*o57*o63-o143*o299*pbx)+o101*o107*o43*o63*(2*o13*(-2*o10*o22*o28*o66*o75+2*o22*o28*o51*o80-2*o10*o51*o66*pax-2*o10*o75*pax+4*o67*o75*pax)+o69*o7*(o308+o313-6*o10*o22*o28*o51+8*o22*o28*o66*o75+8*o51*o66*pax-6*o67*pbx)+2*(o313+o316-2*o10*o22*o28*o51+2*o22*o28*o66*o75+2*o51*o66*pax-2*o67*pbx))+o142*o146*o28*o7*(o142*o192*o37*o43*o57*o63-4*o178*o280*pbx-o143*o192*o299*pbx+o143*o63*(o308+o329-6*o10*o37*o43*o57+8*o151*o37*o43*o75+8*o151*o57*pax-6*o174*pbx)+2*(o316+o329+2*o184*o37*o43*o57+2*o151*o37*o43*o75+2*o151*o57*pax-2*o174*pbx)+2*o32*(-2*o151*o31*o37*o43*o75+4*o37*o43*o57*o76-2*o151*o31*o57*pax+4*o140*o75*pax-2*o31*o75*pax+4*o174*o31*pbx-4*o151*o57*o75*pbx-2*o76*pbx)))

	du[6] = 0.25*G*(o43*(o373+o151*o40*o43)+o28*(o373+o25*o28*o66))+o63*pby-0.5*G*(o20*o28*o366*o9+o20*o366*o43*o9+o28*o34*o63*o9*pby+o34*o43*o63*o9*pby)-0.25*G*(-(o146*o194*o28*o358*o40*o43*o57*o7)-o100*o101*o107*o299*o43*pby-2*o142*o194*o28*o356*o7*(o142*o40*o43*o57*o63-o143*o299*pby)+o101*o107*o43*o63*(2*o13*(-2*o10*o25*o28*o66*o75+2*o25*o28*o51*o80-2*o10*o51*o66*pay-2*o10*o75*pay+4*o67*o75*pay)+o69*o7*(o390+o395-6*o10*o25*o28*o51+8*o25*o28*o66*o75+8*o51*o66*pay-6*o67*pby)+2*(o395+o398-2*o10*o25*o28*o51+2*o25*o28*o66*o75+2*o51*o66*pay-2*o67*pby))+o142*o146*o28*o7*(o142*o192*o40*o43*o57*o63-4*o178*o280*pby-o143*o192*o299*pby+o143*o63*(o390+o411-6*o10*o40*o43*o57+8*o151*o40*o43*o75+8*o151*o57*pay-6*o174*pby)+2*(o398+o411+2*o184*o40*o43*o57+2*o151*o40*o43*o75+2*o151*o57*pay-2*o174*pby)+2*o32*(-2*o151*o31*o40*o43*o75+4*o40*o43*o57*o76-2*o151*o31*o57*pay+4*o140*o75*pay-2*o31*o75*pay+4*o174*o31*pby-4*o151*o57*o75*pby-2*o76*pby)))

	du[7] = -0.25*G*(o22*o443*o451-o37*o445*o478+o43*(o151*o646+o57*o651)+o28*(o51*o639+o634*o66))+0.5*G*(o20*o22*o34*o443*o9-o20*o34*o37*o445*o9)+0.25*G*(-(o100*o101*o107*o37*o445*o63)-o100*o107*o138*o43*o63*o639*o66+o142*o146*o194*o22*o443*o7-o146*o194*o28*o358*o57*o646*o7-2*o194*o28*o356*o532*o57*o63*o646*o7-2*o100*o43*o507*o63*o639*o66*o7*o72+o142*o146*o28*o7*(o142*o192*o57*o63*o646+2*(2*o184*o57*o646-2*o151*o31*o651+o694+o695+2*o151*o646*o75+2*o57*o651*o75)+o143*o63*(-6*o10*o57*o646-6*o151*o31*o651+o694+o695+8*o151*o646*o75+8*o57*o651*o75)+2*o32*(2*o148*o151*o651-2*o151*o31*o646*o75-2*o31*o57*o651*o75+4*o57*o646*o76))+o101*o107*o43*o63*(o69*o7*(-6*o10*o51*o634-6*o31*o639*o66+o669+o670+8*o51*o639*o75+8*o634*o66*o75)+2*o13*(-2*o10*o51*o639*o75-2*o10*o634*o66*o75+4*o639*o66*o76+2*o51*o634*o80)+o101*o639*o66*o7*o90+2*(-2*o10*o51*o634+o669+o670+2*o51*o639*o75+2*o634*o66*o75+2*o639*o66*o94)))

	du[8] = -0.25*G*(o25*o443*o451-o40*o445*o478+o28*(o51*o720+o66*o725)+o43*(o151*o733+o57*o738))+0.5*G*(o20*o25*o34*o443*o9-o20*o34*o40*o445*o9)+0.25*G*(-(o100*o101*o107*o40*o445*o63)+o142*o146*o194*o25*o443*o7-o100*o107*o138*o43*o63*o66*o720-2*o100*o43*o507*o63*o66*o7*o72*o720-o146*o194*o28*o358*o57*o7*o733-2*o194*o28*o356*o532*o57*o63*o7*o733+o142*o146*o28*o7*(o142*o192*o57*o63*o733+2*o32*(2*o148*o151*o738-2*o151*o31*o733*o75-2*o31*o57*o738*o75+4*o57*o733*o76)+2*(2*o184*o57*o733-2*o151*o31*o738+2*o151*o733*o75+2*o57*o738*o75+o781+o782)+o143*o63*(-6*o10*o57*o733-6*o151*o31*o738+8*o151*o733*o75+8*o57*o738*o75+o781+o782))+o101*o107*o43*o63*(o69*o7*(-6*o31*o66*o720-6*o10*o51*o725+8*o51*o720*o75+8*o66*o725*o75+o754+o757)+2*o13*(-2*o10*o51*o720*o75-2*o10*o66*o725*o75+4*o66*o720*o76+2*o51*o725*o80)+o101*o66*o7*o720*o90+2*(-2*o10*o51*o725+2*o51*o720*o75+2*o66*o725*o75+o754+o757+2*o66*o720*o94)))

	u[9] = o20+0.25*G*(o28*o451+o43*o478)-0.25*G*(o100*o101*o107*o43*o63+o142*o146*o194*o28*o7)+o9-0.5*G*(o20*o28*o34*o9+o20*o34*o43*o9)
end

julia_main()

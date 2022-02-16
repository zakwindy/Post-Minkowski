#----------------------------------------------------------------------------
#               Initial condition generator
# File name:    2bodyIC.py
# Creator:      Zackary Windham
# Description:  A code to generate initial conditions for a simple binary syste
#               in either a Newtonian or a post-Minkowskian system.
#----------------------------------------------------------------------------
import argparse
import sys
from math import *
from numpy import *
from scipy import optimize
import pandas

#----------------------------------------------------------------------
#  Set up commandline arguments
#----------------------------------------------------------------------
parser = argparse.ArgumentParser()

parser.add_argument('-o','--order', type=str, default="1",
     help='Post-Minkowskian order of the equations',
     choices=["0","1"])
parser.add_argument('m1', type=float, default=1.0,
       help='Mass of body 1 (Default 1.0)', metavar='NUM')
parser.add_argument('m2', type=float, default=1.0,
       help='Mass of body 2 (Default 1.0)', metavar='NUM')
parser.add_argument('-a', '--semimajor_axis', type=float, default=20.0,
     help='semimajor axis of the binary star system in kilometers(default 20.0)', metavar = 'NUM')
parser.add_argument('-e', '--eccentricity', type=float, default=0.0,
     help='eccentricity of the binary star system', metavar = 'NUM')
parser.add_argument('-w', '--argument_of_periapsis', type=float, default=0.0,
     help='measured from the ascending node to the periapsis in degrees', metavar = 'NUM')
parser.add_argument('-Omega', '--longitude_of_ascending_node', type=float, default=0.0,
     help='measured in degrees, defines where orbital plane passes upward through the reference plane', metavar = 'NUM')
parser.add_argument('-M', '--mean_anomaly', type=float, default=20.0,
     help='measured in degrees, can be used to calculate true anomaly', metavar = 'NUM')
parser.add_argument('-n','--run_number', type=int, default=0,metavar='INT',help="numerical tag for output files")

args = parser.parse_args()

#----------------------------------------------------------------------
#   Set constants and read in parameters
#----------------------------------------------------------------------

C = 1.0            # set relativistic C

C_CGS = 2.998e10
G_CGS = 6.674e-8
mSun_CGS = 1.989e33
AU_CGS = 1.496e14 # 1 AU in cm
KM_CGS = 1e6        # 1 km in cm

# Defined using C = G = 1
G = 1.0
M = mSun_CGS       #units of mass
L = M * (G_CGS / G) * ((C / C_CGS)**2)        #units of length
T = L * C / C_CGS        #units of time
    
m1 = args.m1
m2 = args.m2
Mtot = m1 + m2
mu = m1 * m2 / Mtot

a_KM = args.semimajor_axis
e = args.eccentricity
w_degree = args.argument_of_periapsis
omega_degree = args.longitude_of_ascending_node
Manom_degree = args.mean_anomaly

a = a_KM * KM_CGS / L
i = 0
w = w_degree*pi/180
omega = omega_degree*pi/180
Manom = Manom_degree*pi/180

T_period = 2*pi*sqrt((a**3)/(G*Mtot))

#----------------------------------------------------------------------
#   Calculate initial conditions in Cartesian coordinates
#----------------------------------------------------------------------

u0 = 10        # initializing u0
tol  = 1      # Initalizing iterative tolerance
tol0 = 1e-10

while abs(tol) > tol0:
    u0 = u0 - (-Manom + u0 - e*sin(u0))/(1 - e*cos(u0))
    tol = -Manom + u0 - e*sin(u0)

f     = atan(sin(u0)*sqrt(1-e**2)/(cos(u0)-e))   #use eccentric anomaly to find the true anomaly

# Solving for Polar Coordinates
r     = a*(1-e**2)/(1-e*cos(f))
phi   = omega + atan(tan(w+f)*cos(i))
theta = acos(sin(w+f)*sin(i))

# Transform to Cartesian coordinates
xi=r*sin(theta)*cos(phi)
yi=r*sin(theta)*sin(phi)
zi=r*cos(theta)

x1 = m2*xi/Mtot
y1 = m2*yi/Mtot
z1 = m2*zi/Mtot

x2 = -m1*xi/Mtot
y2 = -m1*yi/Mtot
z2 = -m1*zi/Mtot

gr = a*(1-e**2)*e*sin(f)/((1+e*cos(f))**2)
gtheta = (-1/sin(theta))*cos(w+f)*sin(i)
gphi = (cos(phi-omega)**2)*cos(i)/(cos(w+f)**2)
fdot = sqrt(G*Mtot*((2/r)-(1/a))*1/(gr**2 + (r*gtheta)**2 +(r*sin(theta)*gphi)**2))

rdot = gr*fdot
thetadot = gtheta*fdot
phidot = gphi*fdot

xdot = rdot*sin(theta)*cos(phi) + r*thetadot*cos(theta)*cos(phi) -phidot*r*sin(phi)*sin(theta)
ydot = rdot*sin(theta)*sin(phi) + r*thetadot*cos(theta)*sin(phi) +r*phidot*sin(theta)*cos(phi)
zdot = rdot*cos(theta) - r*thetadot*sin(theta)

#----------------------------------------------------------------------
#  Routines to calculate circular binary initial data
#----------------------------------------------------------------------

#----------------- 0th order-------------------------------------------
def circular_orbit_pm0(r, mu, M):
    px = mu*xdot
    py = mu*ydot
    pz = mu*zdot
    
    p_initial = sqrt(px**2 + py**2 + pz**2)
    return p_initial*r

#----------------- 1st order-------------------------------------------
def f(x):
    o4=1.+e;
    o5=e*e;
    o6=-o5;
    o7=1.+o6;
    o14=x*x;
    o16=pow(a,-3.);
    o17=pow(o4,3.);
    o18=pow(o7,-3.);
    o20=pow(a,-2.);
    o21=o4*o4;
    o22=pow(o7,-2.);
    o23=o14*o20*o21*o22;
    o9=a*a;
    o10=pow(o4,-2.);
    o11=o7*o7;
    o12=o10*o11*o9;
    o19=m1*m1;
    o24=o19+o23;
    o27=m2*m2;
    o28=o23+o27;
    o34=pow(a,-5.);
    o35=pow(o4,5.);
    o36=pow(o7,-5.);
    o37=pow(x,4.);
    o31=1./sqrt(o12);
    o32=sqrt(o24);
    o29=1./sqrt(o28);
    o40=1/o24;
    o44=1/o28;
    o25=1./sqrt(o24);
    o33=sqrt(o28);
    o48=o14*o20*o21*o22*o40;
    o49=o14*o20*o21*o22*o44;
    o50=1.+o48+o49;
    o13=pow(o12,-1.5);
    o53=1/o4;
    o56=sqrt(o23);
    o57=1./sqrt(o23);
    o70=pow(a,-7.);
    o71=pow(o4,7.);
    o72=pow(o7,-7.);
    o73=pow(x,6.);
    o74=pow(o28,-1.5);
    o63=-(o14*o20*o21*o22*o31*o57);
    o64=-(o10*o11*o13*o56*o9);
    o65=o31*o56;
    o66=o63+o64+o65;
    o76=pow(o24,-1.5);
    o81=pow(a,-4.);
    o82=pow(o4,4.);
    o83=pow(o7,-4.);
    o3=1/a;
    o8=1/o7;
    o58=o14*o20*o21*o22*o31*o57;
    o59=o10*o11*o13*o56*o9;
    o60=-(o31*o56);
    o61=o58+o59+o60;
    o92=-8.*o34*o35*o36*o37;
    o94=-4.*o14*o3*o31*o4*o56*o61*o8;
    o95=4.*o14*o3*o31*o4*o56*o66*o8;
    o97=-4.*o34*o35*o36*o37;
    o98=-2.*o14*o3*o31*o4*o56*o61*o8;
    o99=2.*o14*o3*o31*o4*o56*o66*o8;
    o100=o97+o98+o99;
    o101=2.*o100;
    o102=2.*o70*o71*o72*o73;
    o103=2.*o16*o17*o18*o31*o37*o56*o66;
    o104=o102+o103;
    return -(o14*o16*o17*o18*o25)-o14*o16*o17*o18*o29+G*o14*o16*o17*o18*o29*o31*o32*o50+G*o14*o16*o17*o18*o25*o31*o33*o50+a*G*o13*o32*o33*o50*o53*o7+0.5*G*o31*(14.*o14*o16*o17*o18+a*o31*o53*o56*o61*o7-a*o31*o53*o56*o66*o7)+4.*G*o13*o14*o3*o4*o8-0.25*G*(-4.*o13*o16*o17*o18*o25*o29*o37-2.*o25*o34*o35*o36*o37*o56*o66*o74+2.*o25*o31*o70*o71*o72*o73*o74-2.*o29*o34*o35*o36*o37*o56*o66*o76+2.*o29*o31*o70*o71*o72*o73*o76-2.*o25*o29*o31*o37*(o14*o16*o17*o18*o40+a*o31*o40*o53*o56*o66*o7)*o81*o82*o83-2.*o25*o29*o31*o37*(o14*o16*o17*o18*o44+a*o31*o44*o53*o56*o66*o7)*o81*o82*o83+0.25*o25*o29*o31*(o101+2.*o104*o40+4.*o16*o17*o18*o31*o37*o40*o56*o66+4.*o40*o70*o71*o72*o73+o92+o94+o95)+0.25*o25*o29*o31*(o101+2.*o104*o44+4.*o16*o17*o18*o31*o37*o44*o56*o66+4.*o44*o70*o71*o72*o73+o92+o94+o95))-G*o31*o32*o33*(-2.*o14*o16*o17*o18*o40-2.*o14*o16*o17*o18*o44+2.*o34*o35*o36*o37*pow(o24,-2.)+2.*o34*o35*o36*o37*pow(o28,-2.))
    
def fprime(x):
    o4=1.+e;
    o5=e*e;
    o6=-o5;
    o7=1.+o6;
    o19=m1*m1;
    o20=pow(a,-2.);
    o21=o4*o4;
    o22=pow(o7,-2.);
    o23=x*x;
    o24=o20*o21*o22*o23;
    o25=o19+o24;
    o15=pow(a,-5.);
    o16=pow(o4,5.);
    o17=pow(o7,-5.);
    o18=pow(x,3.);
    o28=pow(a,-3.);
    o29=pow(o4,3.);
    o30=pow(o7,-3.);
    o33=m2*m2;
    o34=o24+o33;
    o9=a*a;
    o10=pow(o4,-2.);
    o11=o7*o7;
    o12=o10*o11*o9;
    o42=pow(a,-7.);
    o43=pow(o4,7.);
    o44=pow(o7,-7.);
    o45=pow(x,5.);
    o39=1./sqrt(o12);
    o40=sqrt(o25);
    o37=1./sqrt(o34);
    o48=pow(o25,-2.);
    o50=1/o25;
    o60=pow(a,-4.);
    o61=pow(o4,4.);
    o62=pow(o7,-4.);
    o54=pow(o34,-2.);
    o56=1/o34;
    o31=1./sqrt(o25);
    o41=sqrt(o34);
    o63=-2.*o18*o48*o60*o61*o62;
    o64=2.*o20*o21*o22*o50*x;
    o65=-2.*o18*o54*o60*o61*o62;
    o66=2.*o20*o21*o22*o56*x;
    o67=o63+o64+o65+o66;
    o13=pow(o12,-1.5);
    o72=pow(x,4.);
    o73=2.*o15*o16*o17*o48*o72;
    o74=-2.*o23*o28*o29*o30*o50;
    o75=2.*o15*o16*o17*o54*o72;
    o76=-2.*o23*o28*o29*o30*o56;
    o77=o73+o74+o75+o76;
    o35=pow(o34,-1.5);
    o80=o20*o21*o22*o23*o50;
    o81=o20*o21*o22*o23*o56;
    o82=1.+o80+o81;
    o3=1/a;
    o8=1/o7;
    o26=pow(o25,-1.5);
    o70=1/o4;
    o94=1./sqrt(o24);
    o91=sqrt(o24);
    o92=pow(o24,-1.5);
    o116=pow(a,-9.);
    o117=pow(o4,9.);
    o118=pow(o7,-9.);
    o119=pow(x,7.);
    o93=o18*o39*o60*o61*o62*o92;
    o95=-(o13*o94*x);
    o96=-(o20*o21*o22*o39*o94*x);
    o97=o93+o95+o96;
    o120=pow(o34,-2.5);
    o109=-(o20*o21*o22*o23*o39*o94);
    o110=-(o10*o11*o13*o9*o91);
    o111=o39*o91;
    o112=o109+o110+o111;
    o125=pow(o25,-2.5);
    o146=pow(a,-6.);
    o147=pow(o4,6.);
    o148=pow(o7,-6.);
    o149=o23*o28*o29*o30*o50;
    o150=a*o112*o39*o50*o7*o70*o91;
    o151=o149+o150;
    o162=o23*o28*o29*o30*o56;
    o163=a*o112*o39*o56*o7*o70*o91;
    o164=o162+o163;
    o99=-(o18*o39*o60*o61*o62*o92);
    o100=o13*o94*x;
    o101=o20*o21*o22*o39*o94*x;
    o102=o100+o101+o99;
    o104=o20*o21*o22*o23*o39*o94;
    o105=o10*o11*o13*o9*o91;
    o106=-(o39*o91);
    o107=o104+o105+o106;
    o196=pow(x,6.);
    o197=2.*o196*o42*o43*o44;
    o198=2.*o112*o28*o29*o30*o39*o72*o91;
    o199=o197+o198;
    o203=-8.*o15*o16*o17*o72;
    o204=4.*o196*o42*o43*o44*o50;
    o205=-4.*o107*o23*o3*o39*o4*o8*o91;
    o206=4.*o112*o23*o3*o39*o4*o8*o91;
    o207=4.*o112*o28*o29*o30*o39*o50*o72*o91;
    o208=-4.*o15*o16*o17*o72;
    o209=-2.*o107*o23*o3*o39*o4*o8*o91;
    o210=2.*o112*o23*o3*o39*o4*o8*o91;
    o211=o208+o209+o210;
    o212=2.*o211;
    o213=2.*o199*o50;
    o214=o203+o204+o205+o206+o207+o212+o213;
    o168=-32.*o15*o16*o17*o18;
    o171=4.*o23*o3*o39*o4*o8*o91*o97;
    o173=-4.*o102*o23*o3*o39*o4*o8*o91;
    o174=-4.*o107*o18*o28*o29*o30*o39*o94;
    o175=-8.*o107*o3*o39*o4*o8*o91*x;
    o176=4.*o112*o18*o28*o29*o30*o39*o94;
    o177=8.*o112*o3*o39*o4*o8*o91*x;
    o181=-16.*o15*o16*o17*o18;
    o182=2.*o23*o3*o39*o4*o8*o91*o97;
    o183=-2.*o102*o23*o3*o39*o4*o8*o91;
    o184=-2.*o107*o18*o28*o29*o30*o39*o94;
    o185=-4.*o107*o3*o39*o4*o8*o91*x;
    o186=2.*o112*o18*o28*o29*o30*o39*o94;
    o187=4.*o112*o3*o39*o4*o8*o91*x;
    o188=o181+o182+o183+o184+o185+o186+o187;
    o189=2.*o188;
    o190=12.*o42*o43*o44*o45;
    o191=2.*o28*o29*o30*o39*o72*o91*o97;
    o192=2.*o112*o15*o16*o17*o39*o45*o94;
    o193=8.*o112*o18*o28*o29*o30*o39*o91;
    o194=o190+o191+o192+o193;
    o227=4.*o196*o42*o43*o44*o56;
    o228=4.*o112*o28*o29*o30*o39*o56*o72*o91;
    o229=2.*o199*o56;
    o230=o203+o205+o206+o212+o227+o228+o229;
    return o15*o16*o17*o18*o26+o15*o16*o17*o18*o35+G*o23*o28*o29*o30*o37*o39*o40*o67+G*o23*o28*o29*o30*o31*o39*o41*o67+a*G*o13*o40*o41*o67*o7*o70+2.*G*o15*o16*o17*o18*o31*o37*o39*o82-G*o15*o16*o17*o18*o35*o39*o40*o82-G*o15*o16*o17*o18*o26*o39*o41*o82-2.*o28*o29*o30*o31*x-2.*o28*o29*o30*o37*x-G*o20*o21*o22*o37*o39*o40*o77*x-G*o20*o21*o22*o31*o39*o41*o77*x+8.*G*o13*o3*o4*o8*x+2.*G*o28*o29*o30*o37*o39*o40*o82*x+2.*G*o28*o29*o30*o31*o39*o41*o82*x+G*o13*o3*o37*o4*o40*o8*o82*x+G*o13*o3*o31*o4*o41*o8*o82*x+0.5*G*o39*(a*o102*o39*o7*o70*o91-a*o39*o7*o70*o91*o97+28.*o28*o29*o30*x+o107*o3*o39*o4*o8*o94*x-o112*o3*o39*o4*o8*o94*x)-0.25*G*(-16.*o13*o18*o28*o29*o30*o31*o37-6.*o116*o117*o118*o119*o120*o31*o39-4.*o116*o117*o118*o119*o26*o35*o39-6.*o116*o117*o118*o119*o125*o37*o39+4.*o13*o15*o16*o17*o31*o35*o45+4.*o13*o15*o16*o17*o26*o37*o45+2.*o146*o147*o148*o151*o31*o35*o39*o45+2.*o146*o147*o148*o164*o31*o35*o39*o45+2.*o146*o147*o148*o151*o26*o37*o39*o45+2.*o146*o147*o148*o164*o26*o37*o39*o45+12.*o31*o35*o39*o42*o43*o44*o45+12.*o26*o37*o39*o42*o43*o44*o45-8.*o151*o18*o31*o37*o39*o60*o61*o62-8.*o164*o18*o31*o37*o39*o60*o61*o62-8.*o112*o15*o16*o17*o18*o31*o35*o91-8.*o112*o15*o16*o17*o18*o26*o37*o91+6.*o112*o120*o31*o42*o43*o44*o45*o91+4.*o112*o26*o35*o42*o43*o44*o45*o91+6.*o112*o125*o37*o42*o43*o44*o45*o91-2.*o112*o31*o35*o42*o43*o44*o45*o94-2.*o112*o26*o37*o42*o43*o44*o45*o94-2.*o15*o16*o17*o31*o35*o72*o91*o97-2.*o15*o16*o17*o26*o37*o72*o91*o97-0.25*o20*o21*o214*o22*o31*o35*o39*x-0.25*o20*o21*o22*o230*o31*o35*o39*x-0.25*o20*o21*o214*o22*o26*o37*o39*x-0.25*o20*o21*o22*o230*o26*o37*o39*x+0.25*o31*o37*o39*(o168+o171+o173+o174+o175+o176+o177+o189-8.*o116*o117*o118*o119*o48+2.*o194*o50+24.*o42*o43*o44*o45*o50-8.*o112*o15*o16*o17*o39*o45*o48*o91+16.*o112*o18*o28*o29*o30*o39*o50*o91+4.*o112*o15*o16*o17*o39*o45*o50*o94+4.*o28*o29*o30*o39*o50*o72*o91*o97-4.*o199*o20*o21*o22*o48*x)+0.25*o31*o37*o39*(o168+o171+o173+o174+o175+o176+o177+o189-8.*o116*o117*o118*o119*o54+2.*o194*o56+24.*o42*o43*o44*o45*o56-8.*o112*o15*o16*o17*o39*o45*o54*o91+16.*o112*o18*o28*o29*o30*o39*o56*o91+4.*o112*o15*o16*o17*o39*o45*o56*o94+4.*o28*o29*o30*o39*o56*o72*o91*o97-4.*o199*o20*o21*o22*o54*x)-2.*o31*o37*o39*o60*o61*o62*o72*(-2.*o15*o16*o17*o18*o48+a*o39*o50*o7*o70*o91*o97+2.*o28*o29*o30*o50*x-2.*o112*o3*o39*o4*o48*o8*o91*x+o112*o3*o39*o4*o50*o8*o94*x)-2.*o31*o37*o39*o60*o61*o62*o72*(-2.*o15*o16*o17*o18*o54+a*o39*o56*o7*o70*o91*o97+2.*o28*o29*o30*o56*x-2.*o112*o3*o39*o4*o54*o8*o91*x+o112*o3*o39*o4*o56*o8*o94*x))-G*o39*o40*o41*(12.*o15*o16*o17*o18*o48+12.*o15*o16*o17*o18*o54-4.*o28*o29*o30*o50*x-4.*o28*o29*o30*o56*x-8.*o42*o43*o44*o45*pow(o25,-3.)-8.*o42*o43*o44*o45*pow(o34,-3.));

    
def fprime2(x):
    o4=1.+e;
    o5=e*e;
    o6=-o5;
    o7=1.+o6;
    o17=pow(a,-2.);
    o18=o4*o4;
    o20=x*x;
    o21=m1*m1;
    o22=pow(o7,-2.);
    o23=o17*o18*o20*o22;
    o24=o21+o23;
    o27=pow(o7,-3.);
    o25=pow(o24,-1.5);
    o15=pow(a,-3.);
    o16=pow(o4,3.);
    o19=pow(o7,-5.);
    o40=m2*m2;
    o41=o23+o40;
    o30=pow(a,-4.);
    o31=pow(o4,4.);
    o32=pow(o7,-4.);
    o42=pow(o41,-1.5);
    o9=a*a;
    o10=pow(o4,-2.);
    o11=o7*o7;
    o12=o10*o11*o9;
    o54=pow(a,-7.);
    o55=pow(o4,7.);
    o56=pow(o7,-7.);
    o57=pow(x,5.);
    o60=pow(a,-5.);
    o61=pow(o4,5.);
    o62=pow(x,3.);
    o44=1./sqrt(o41);
    o28=1./sqrt(o24);
    o63=pow(o24,-2.);
    o65=1/o24;
    o80=pow(x,4.);
    o69=pow(o41,-2.);
    o71=1/o41;
    o76=sqrt(o41);
    o74=sqrt(o24);
    o58=pow(o24,-3.);
    o67=pow(o41,-3.);
    o53=1./sqrt(o12);
    o13=pow(o12,-1.5);
    o137=1./sqrt(o23);
    o138=pow(o23,-1.5);
    o148=sqrt(o23);
    o144=-(o138*o20*o30*o31*o32);
    o145=o137*o17*o18*o22;
    o146=o144+o145;
    o136=1/o4;
    o155=-4.*o138*o17*o18*o20*o22*o53;
    o156=2.*o137*o53;
    o157=pow(o23,-2.5);
    o158=3.*o157*o20*o30*o31*o32;
    o159=-(o138*o17*o18*o22);
    o160=o158+o159;
    o161=o160*o20*o53;
    o162=o155+o156+o161;
    o87=-(o20*o25*o30*o31*o32);
    o88=o17*o18*o22*o28;
    o89=o87+o88;
    o187=o17*o18*o20*o22*o65;
    o188=o17*o18*o20*o22*o71;
    o189=1.+o187+o188;
    o198=-2.*o30*o31*o32*o62*o63;
    o199=2.*o17*o18*o22*o65*x;
    o200=-2.*o30*o31*o32*o62*o69;
    o201=2.*o17*o18*o22*o71*x;
    o202=o198+o199+o200+o201;
    o46=pow(o41,-2.5);
    o47=3.*o20*o30*o31*o32*o46;
    o48=-(o17*o18*o22*o42);
    o49=o47+o48;
    o108=8.*o20*o30*o31*o32*o58;
    o109=-2.*o17*o18*o22*o63;
    o110=o108+o109;
    o125=8.*o20*o30*o31*o32*o67;
    o126=-2.*o17*o18*o22*o69;
    o127=o125+o126;
    o33=pow(o24,-2.5);
    o34=3.*o20*o30*o31*o32*o33;
    o35=-(o17*o18*o22*o25);
    o36=o34+o35;
    o91=-(o20*o30*o31*o32*o42);
    o92=o17*o18*o22*o44;
    o93=o91+o92;
    o209=-8.*o17*o18*o20*o22*o63;
    o210=2.*o65;
    o211=o110*o20;
    o212=o209+o210+o211;
    o213=o17*o18*o212*o22;
    o214=-8.*o17*o18*o20*o22*o69;
    o215=2.*o71;
    o216=o127*o20;
    o217=o214+o215+o216;
    o218=o17*o18*o217*o22;
    o219=o213+o218;
    o233=o202*o76;
    o234=o17*o18*o189*o22*o44*x;
    o235=o233+o234;
    o237=2.*o17*o18*o202*o22*o44*x;
    o238=o189*o93;
    o239=o219*o76;
    o240=o237+o238+o239;
    o147=-(o137*o17*o18*o20*o22*o53);
    o149=-(o10*o11*o13*o148*o9);
    o150=o148*o53;
    o151=o147+o149+o150;
    o139=o138*o30*o31*o32*o53*o62;
    o140=-(o13*o137*x);
    o141=-(o137*o17*o18*o22*o53*x);
    o142=o139+o140+o141;
    o253=-(o148*o151*o17*o18*o22*o25*o42*x);
    o3=1/a;
    o8=1/o7;
    o173=o137*o17*o18*o20*o22*o53;
    o174=o10*o11*o13*o148*o9;
    o175=-(o148*o53);
    o176=o173+o174+o175;
    o263=pow(x,6.);
    o262=-8.*o19*o60*o61*o80;
    o265=-4.*o148*o176*o20*o3*o4*o53*o8;
    o266=4.*o148*o151*o20*o3*o4*o53*o8;
    o268=-4.*o19*o60*o61*o80;
    o269=-2.*o148*o176*o20*o3*o4*o53*o8;
    o270=2.*o148*o151*o20*o3*o4*o53*o8;
    o271=o268+o269+o270;
    o272=2.*o271;
    o273=2.*o263*o54*o55*o56;
    o274=2.*o148*o15*o151*o16*o27*o53*o80;
    o275=o273+o274;
    o286=-(o17*o18*o22*o25*o42*x);
    o168=-(o138*o30*o31*o32*o53*o62);
    o169=o13*o137*x;
    o170=o137*o17*o18*o22*o53*x;
    o171=o168+o169+o170;
    o264=4.*o263*o54*o55*o56*o65;
    o267=4.*o148*o15*o151*o16*o27*o53*o65*o80;
    o276=2.*o275*o65;
    o277=o262+o264+o265+o266+o267+o272+o276;
    o327=-32.*o19*o60*o61*o62;
    o328=pow(a,-9.);
    o329=pow(o4,9.);
    o330=pow(o7,-9.);
    o331=pow(x,7.);
    o334=4.*o142*o148*o20*o3*o4*o53*o8;
    o336=-4.*o148*o171*o20*o3*o4*o53*o8;
    o337=-4.*o137*o15*o16*o176*o27*o53*o62;
    o338=-8.*o148*o176*o3*o4*o53*o8*x;
    o339=4.*o137*o15*o151*o16*o27*o53*o62;
    o340=8.*o148*o151*o3*o4*o53*o8*x;
    o344=-16.*o19*o60*o61*o62;
    o345=2.*o142*o148*o20*o3*o4*o53*o8;
    o346=-2.*o148*o171*o20*o3*o4*o53*o8;
    o347=-2.*o137*o15*o16*o176*o27*o53*o62;
    o348=-4.*o148*o176*o3*o4*o53*o8*x;
    o349=2.*o137*o15*o151*o16*o27*o53*o62;
    o350=4.*o148*o151*o3*o4*o53*o8*x;
    o351=o344+o345+o346+o347+o348+o349+o350;
    o352=2.*o351;
    o353=12.*o54*o55*o56*o57;
    o354=2.*o142*o148*o15*o16*o27*o53*o80;
    o355=2.*o137*o151*o19*o53*o57*o60*o61;
    o356=8.*o148*o15*o151*o16*o27*o53*o62;
    o357=o353+o354+o355+o356;
    o279=4.*o263*o54*o55*o56*o71;
    o280=4.*o148*o15*o151*o16*o27*o53*o71*o80;
    o281=2.*o275*o71;
    o282=o262+o265+o266+o272+o279+o280+o281;
    o291=pow(o41,-3.5);
    o292=15.*o20*o291*o30*o31*o32;
    o293=-3.*o17*o18*o22*o46;
    o294=o292+o293;
    o153=-(o10*o11*o13*o146*o9);
    o154=o146*o53;
    o163=-(o162*o17*o18*o22);
    o164=o153+o154+o163;
    o305=pow(o24,-3.5);
    o306=15.*o20*o30*o305*o31*o32;
    o307=-3.*o17*o18*o22*o33;
    o308=o306+o307;
    o416=o15*o16*o20*o27*o65;
    o417=a*o136*o148*o151*o53*o65*o7;
    o418=o416+o417;
    o427=-2.*o19*o60*o61*o62*o63;
    o428=2.*o15*o16*o27*o65*x;
    o429=a*o136*o142*o148*o53*o65*o7;
    o430=-2.*o148*o151*o3*o4*o53*o63*o8*x;
    o431=o137*o151*o3*o4*o53*o65*o8*x;
    o432=o427+o428+o429+o430+o431;
    o419=-8.*o17*o18*o22*o25*o80;
    o420=12.*o20*o28;
    o421=o36*o80;
    o422=o419+o420+o421;
    o424=-(o17*o18*o22*o25*o57);
    o425=4.*o28*o62;
    o426=o424+o425;
    o458=o15*o16*o20*o27*o71;
    o459=a*o136*o148*o151*o53*o7*o71;
    o460=o458+o459;
    o462=-2.*o19*o60*o61*o62*o69;
    o463=2.*o15*o16*o27*o71*x;
    o464=a*o136*o142*o148*o53*o7*o71;
    o465=-2.*o148*o151*o3*o4*o53*o69*o8*x;
    o466=o137*o151*o3*o4*o53*o71*o8*x;
    o467=o462+o463+o464+o465+o466;
    o332=-8.*o328*o329*o330*o331*o63;
    o333=24.*o54*o55*o56*o57*o65;
    o335=4.*o142*o148*o15*o16*o27*o53*o65*o80;
    o341=-8.*o148*o151*o19*o53*o57*o60*o61*o63;
    o342=4.*o137*o151*o19*o53*o57*o60*o61*o65;
    o343=16.*o148*o15*o151*o16*o27*o53*o62*o65;
    o358=2.*o357*o65;
    o359=-4.*o17*o18*o22*o275*o63*x;
    o360=o327+o332+o333+o334+o335+o336+o337+o338+o339+o340+o341+o342+o343+o352+o358+o359;
    o178=o10*o11*o13*o146*o9;
    o179=-(o146*o53);
    o180=o162*o17*o18*o22;
    o181=o178+o179+o180;
    o504=o142*o148;
    o505=o137*o151*o17*o18*o22*x;
    o506=o504+o505;
    o508=2.*o137*o142*o17*o18*o22*x;
    o509=o146*o151;
    o510=o148*o164;
    o511=o508+o509+o510;
    o503=2.*o148*o151*o53;
    o507=4.*o506*o53*x;
    o512=o20*o511*o53;
    o513=o503+o507+o512;
    o534=2.*o148*o176*o53;
    o535=o148*o171;
    o536=o137*o17*o176*o18*o22*x;
    o537=o535+o536;
    o538=4.*o53*o537*x;
    o539=2.*o137*o17*o171*o18*o22*x;
    o540=o146*o176;
    o541=o148*o181;
    o542=o539+o540+o541;
    o543=o20*o53*o542;
    o544=o534+o538+o543;
    o365=-8.*o328*o329*o330*o331*o69;
    o366=24.*o54*o55*o56*o57*o71;
    o367=4.*o142*o148*o15*o16*o27*o53*o71*o80;
    o368=-8.*o148*o151*o19*o53*o57*o60*o61*o69;
    o369=4.*o137*o151*o19*o53*o57*o60*o61*o71;
    o370=16.*o148*o15*o151*o16*o27*o53*o62*o71;
    o371=2.*o357*o71;
    o372=-4.*o17*o18*o22*o275*o69*x;
    o373=o327+o334+o336+o337+o338+o339+o340+o352+o365+o366+o367+o368+o369+o370+o371+o372;
    o495=-96.*o19*o20*o60*o61;
    o514=4.*o3*o4*o513*o8;
    o515=8.*o137*o17*o18*o22*o80;
    o516=12.*o148*o20;
    o517=o146*o80;
    o518=o515+o516+o517;
    o520=o137*o17*o18*o22*o57;
    o521=4.*o148*o62;
    o522=o520+o521;
    o545=-4.*o3*o4*o544*o8;
    o546=60.*o54*o55*o56*o80;
    o547=12.*o148*o151*o20*o53;
    o548=8.*o506*o53*o62;
    o549=o511*o53*o80;
    o550=o547+o548+o549;
    o551=2.*o15*o16*o27*o550;
    o552=o546+o551;
    o554=-48.*o19*o20*o60*o61;
    o555=2.*o3*o4*o513*o8;
    o556=-2.*o3*o4*o544*o8;
    o557=o554+o555+o556;
    o558=2.*o557;
    return -(o15*o16*(-4.*o17*o18*o19*o20*o25+2.*o27*o28+o20*o27*o36))-o15*o16*(-4.*o17*o18*o19*o20*o42+2.*o27*o44+o20*o27*o49)+8.*G*o13*o3*o4*o8+a*o13*o136*o7*(G*o240*o74+G*o189*o76*o89+2.*G*o17*o18*o22*o235*o28*x)+G*o15*o16*o27*o53*(o20*o240*o28+o189*(-4.*o17*o18*o20*o22*o25+2.*o28+o20*o36)*o76+2.*o235*(-(o17*o18*o22*o25*o62)+2.*o28*x))+0.5*G*o53*(28.*o15*o16*o27-a*o136*o7*(o146*o151*o53+o148*o164*o53+2.*o137*o142*o17*o18*o22*o53*x)+a*o136*o7*(o146*o176*o53+o148*o181*o53+2.*o137*o17*o171*o18*o22*o53*x))+G*o15*o16*o27*o53*(o189*o44*(4.*o17*o18*o20*o22*o28+2.*o74+o20*o89)+o20*o74*(o219*o44+o189*o49-2.*o17*o18*o202*o22*o42*x)+2.*(o202*o44-o17*o18*o189*o22*o42*x)*(o17*o18*o22*o28*o62+2.*o74*x))-0.25*G*(0.25*o277*o36*o44*o53+0.25*o282*o36*o44*o53-24.*o148*o151*o19*o20*o28*o42*o60*o61-24.*o148*o151*o19*o20*o25*o44*o60*o61-0.5*o17*o18*o22*o25*o53*x*(o360*o44-o17*o18*o22*o277*o42*x)-0.5*o17*o18*o22*o25*o53*x*(o373*o44-o17*o18*o22*o282*o42*x)-16.*o19*o60*o61*o62*(o253+o142*o148*o25*o44+o137*o151*o17*o18*o22*o25*o44*x-3.*o148*o151*o17*o18*o22*o33*o44*x)-16.*o19*o60*o61*o62*(o253+o142*o148*o28*o42+o137*o151*o17*o18*o22*o28*o42*x-3.*o148*o151*o17*o18*o22*o28*o46*x)-2.*o19*o60*o61*o80*(o151*(o146*o25+o148*o308-6.*o137*o20*o30*o31*o32*o33)*o44+o148*o25*(o164*o44+o151*o49-2.*o142*o17*o18*o22*o42*x)+2.*(o137*o17*o18*o22*o25*x-3.*o148*o17*o18*o22*o33*x)*(o142*o44-o151*o17*o18*o22*o42*x))-4.*o15*o16*o27*(12.*o13*o20*o28*o44+o13*(2.*o20*o25*o30*o31*o32*o42+o36*o44+o28*o49)*o80+8.*o13*o62*(-(o17*o18*o22*o28*o42*x)-o17*o18*o22*o25*o44*x))+2.*o54*o55*o56*(o263*(6.*o20*o30*o31*o32*o33*o42+o308*o44+o25*o49)*o53+30.*o25*o44*o53*o80+12.*o53*o57*(o286-3.*o17*o18*o22*o33*o44*x))-2.*o19*o60*o61*o80*(o151*(o146*o28-2.*o137*o20*o25*o30*o31*o32+o148*o36)*o42+o148*o28*(o151*o294+o164*o42-6.*o142*o17*o18*o22*o46*x)+2.*(-(o148*o17*o18*o22*o25*x)+o137*o17*o18*o22*o28*x)*(o142*o42-3.*o151*o17*o18*o22*o46*x))+2.*o54*o55*o56*(o263*(o28*o294+o36*o42+6.*o20*o25*o30*o31*o32*o46)*o53+30.*o28*o42*o53*o80+12.*o53*o57*(o286-3.*o17*o18*o22*o28*o46*x))+0.25*o28*o53*(o277*o49-2.*o17*o18*o22*o360*o42*x+o44*(2.*o110*o275+o495+o514+o545+o558+2.*o552*o65+4.*o54*o55*(o110*o263*o56-24.*o17*o18*o263*o330*o63+30.*o56*o65*o80)-8.*o17*o18*o22*o357*o63*x+4.*o15*o16*o27*o53*(o151*o518*o65+o148*o80*(o110*o151+o164*o65-4.*o142*o17*o18*o22*o63*x)+2.*o522*(o142*o65-2.*o151*o17*o18*o22*o63*x))))+0.25*o28*o53*(o282*o49-2.*o17*o18*o22*o373*o42*x+o44*(2.*o127*o275+o495+o514+o545+o558+2.*o552*o71+4.*o54*o55*(o127*o263*o56-24.*o17*o18*o263*o330*o69+30.*o56*o71*o80)-8.*o17*o18*o22*o357*o69*x+4.*o15*o16*o27*o53*(o151*o518*o71+o148*o80*(o127*o151+o164*o71-4.*o142*o17*o18*o22*o69*x)+2.*o522*(o142*o71-2.*o151*o17*o18*o22*o69*x))))-2.*o30*o31*o32*o53*(o418*o422*o44+2.*o426*(o432*o44-o17*o18*o22*o418*o42*x)+o28*o80*(o418*o49-2.*o17*o18*o22*o42*o432*x+o44*(o15*o16*o212*o27+a*o136*o53*o7*(o148*o164*o65+o151*(o110*o148-4.*o137*o20*o30*o31*o32*o63+o146*o65)+2.*o142*(-2.*o148*o17*o18*o22*o63*x+o137*o17*o18*o22*o65*x)))))-2.*o30*o31*o32*o53*(o422*o44*o460+2.*o426*(o44*o467-o17*o18*o22*o42*o460*x)+o28*o80*(o460*o49-2.*o17*o18*o22*o42*o467*x+o44*(o15*o16*o217*o27+a*o136*o53*o7*(o148*o164*o71+o151*(o127*o148-4.*o137*o20*o30*o31*o32*o69+o146*o71)+2.*o142*(-2.*o148*o17*o18*o22*o69*x+o137*o17*o18*o22*o71*x))))))-G*o53*((-2.*o15*o16*o20*o27*o65-2.*o15*o16*o20*o27*o71+2.*o19*o60*o61*o63*o80+2.*o19*o60*o61*o69*o80)*(2.*o20*o28*o30*o31*o32*o44+o76*o89+o74*o93)+2.*(-8.*o54*o55*o56*o57*o58+12.*o19*o60*o61*o62*o63-8.*o54*o55*o56*o57*o67+12.*o19*o60*o61*o62*o69-4.*o15*o16*o27*o65*x-4.*o15*o16*o27*o71*x)*(o17*o18*o22*o44*o74*x+o17*o18*o22*o28*o76*x)+o74*o76*(-2.*o15*o16*(o110*o20*o27-8.*o17*o18*o19*o20*o63+2.*o27*o65)-2.*o15*o16*(o127*o20*o27-8.*o17*o18*o19*o20*o69+2.*o27*o71)+2.*o60*o61*(12.*o19*o20*o63-32.*o17*o18*o56*o58*o80+o19*o80*(-4.*o17*o18*o22*o58+24.*o20*o30*o31*o32*pow(o24,-4.)))+2.*o60*o61*(12.*o19*o20*o69-32.*o17*o18*o56*o67*o80+o19*o80*(-4.*o17*o18*o22*o67+24.*o20*o30*o31*o32*pow(o41,-4.)))));

    
def circular_orbit_pm1(r, mu, M):
    sol = optimize.root_scalar(f, x0=circular_orbit_pm0(r, mu, M), fprime=fprime, fprime2=fprime2, method='newton')
    return sol.root

if args.order == "0":
    pphi = circular_orbit_pm0(r, mu, M)
    order = 0
elif args.order == "1":
    pphi = circular_orbit_pm1(r, mu, M)
    order = 1
else:
    print ("Unknown equation order = " + args.order)
    sys.exit(2)
    
px1 = -pphi/(r)*sin(phi)*sin(theta)
py1 = pphi/(r)*cos(phi)*sin(theta)
pz1 = pphi/(r)*cos(theta)
px2 = pphi/(r)*sin(phi)*sin(theta)
py2 = -pphi/(r)*cos(phi)*sin(theta)
pz2 = -pphi/(r)*cos(theta)

CX = (m1*x1+m2*x2)/Mtot
CY = (m1*y1+m2*y2)/Mtot
CZ = (m1*z1+m2*z2)/Mtot

CPX = (px1+px2)/Mtot
CPY = (py1+py2)/Mtot
CPZ = (pz1+pz2)/Mtot

#position transformation
x1 = x1-CX
y1 = y1-CY
z1 = z1-CZ
x2 = x2-CX
y2 = y2-CY
z2 = z2-CZ

#momentum transformation
px1 = px1-CPX
py1 = py1-CPY
pz1 = pz1-CPZ
px2 = px2-CPX
py2 = py2-CPY
pz2 = pz2-CPZ

filename = '2_body_ICfile' + str(args.run_number)
f = open(filename,'w')

s0 = str(G) + " " + str(M) + " " + str(L) + " " + str(T) + " " + str(T_period) + " 0 0\n"
s1 = str(m1) + " " + str(x1) + " " + str(y1) + " " + str(z1) + " " + str(px1) + " " +  str(py1) + " " + str(pz1) + "\n"
s2 = str(m2) + " " + str(x2) + " " + str(y2) + " " + str(z2) + " " + str(px2) + " " +  str(py2) + " " + str(pz2) + "\n"

f.write(s0)
f.write(s1)
f.write(s2)

print("period = ", str(T_period*T/(365*24*3600)), " years.")

f.close()

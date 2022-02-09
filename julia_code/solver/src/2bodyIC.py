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
M = m1 + m2
mu = m1 * m2 / M

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

#----------------------------------------------------------------------
#   Calculate initial conditions in Cartesian coordinates
#----------------------------------------------------------------------

u0 = 10        # initializing u0
tol  = 1      # Initalizing iterative tolerance
tol0 = 1e-10

while abs(tol) > tol0:
    u0 = u0 - (-M + u0 - e*sin(u0))/(1 - e*cos(u0))
    tol = -M + u0 - e*sin(u0)

f     = atan(sin(u0)*sqrt(1-e**2)/(cos(u0)-e))   #use eccentric anomaly to find the true anomaly

# Solving for Polar Coordinates
r     = a*(1-e**2)/(1-e*cos(f))
phi   = omega + atan(tan(w+f)*cos(i))
theta = acos(sin(w+f)*sin(i))

# Transform to Cartesian coordinates
xi=r*sin(theta)*cos(phi)
yi=r*sin(theta)*sin(phi)
zi=r*cos(theta)

x1 = m2*xi/M
y1 = m2*yi/M
z1 = m2*zi/M

x2 = -m1*xi/M
y2 = -m1*yi/M
z2 = -m1*zi/M

#----------------------------------------------------------------------
#  Routines to calculate circular binary initial data
#----------------------------------------------------------------------

#----------------- 0th order-------------------------------------------
def circular_orbit_pm0(r, mu, M):
    q  = G * m1 * m2 * r * mu
    return sqrt(q)

#----------------- 1st order-------------------------------------------
def f(x):
    o29=r*r
    o46=x*x
    o44=m1*m1
    o53=m2*m2
    o60=(r**-2)
    o61=o46*o60
    o45=o29*o44
    o47=o45+o46
    o54=o29*o53
    o58=o46+o54
    o62=o44+o61
    o104=o53+o61
    o108=sqrt(o62)
    o109=sqrt(o104)
    o113=(m1**4)
    o132=(r**4)
    o114=(m2**4)
    
    return (1*(r**-5)*(-((o108+o109)*o46*o47*o58*(o29**1.5))+G*(6*o44*o46*o53*(2*o108*o109+o44+o53)*(r**6)+o113*o114*(r**8)+(4*o113*o132+4*o114*o132+12*o108*o109*o132*o53+o44*(12*o108*o109*o132+27*o132*o53))*(x**4)+(12*o108*o109*o29+18*o29*o44+18*o29*o53)*(x**6)+12*(x**8))))/(o47*o58*sqrt(o104)*sqrt(o29)*sqrt(o62))
    
def fprime(x):
    o29=r*r
    o46=x*x
    o44=m1*m1
    o53=m2*m2
    o60=(r**-2)
    o61=o46*o60
    o45=o29*o44
    o47=o45+o46
    o54=o29*o53
    o58=o46+o54
    o62=o44+o61
    o103=1/sqrt(o62)
    o104=o53+o61
    o105=1/sqrt(o104)
    o106=(o29**1.5)
    o112=(x**3)
    o113=sqrt(o62)
    o114=sqrt(o104)
    o115=o113+o114
    o129=(r**6)
    o165=(r**4)
    o35=1/sqrt(o29)
    o52=1/o47
    o59=1/o58
    o166=(m1**4)
    o168=(m2**4)
    o134=2*o113*o114
    o135=o134+o44+o53
    o122=(x**6)
    o139=18*o29*o44
    o147=18*o29*o53
    o152=12*o113*o114*o29
    o153=o139+o147+o152
    o155=(x**4)
    o167=4*o165*o166
    o169=4*o165*o168
    o170=12*o113*o114*o165*o53
    o171=27*o165*o53
    o172=12*o113*o114*o165
    o195=o171+o172
    o196=o195*o44
    o197=o167+o169+o170+o196
    o206=(r**-7)
    o208=-(o106*o115*o46*o47*o58)
    o210=(r**8)
    o211=o166*o168*o210
    o212=(x**8)
    o213=12*o212
    o214=6*o129*o135*o44*o46*o53
    o215=o122*o153
    o216=o155*o197
    o217=o211+o213+o214+o215+o216
    o218=G*o217
    o219=o208+o218
    o8=(r**-5)
    
    return -(o103*o206*o219*o35*o52*o59*x*(o104**-1.5))-2*o103*o105*o219*o35*o59*o8*x*(o47**-2)-2*o103*o105*o219*o35*o52*o8*x*(o58**-2)-o105*o206*o219*o35*o52*o59*x*(o62**-1.5)+o103*o105*o35*o52*o59*o8*(-2*o106*o112*o115*o47-2*o106*o112*o115*o58-2*o106*o115*o47*o58*x-o106*o46*o47*o58*(o103*o60*x+o105*o60*x)+G*(4*o112*o197+12*o129*o135*o44*o53*x+o122*(12*o105*o113*x+12*o103*o114*x)+6*o129*o44*o46*o53*(2*o105*o113*o60*x+2*o103*o114*o60*x)+o155*(12*o105*o113*o29*o53*x+12*o103*o114*o29*o53*x+o44*(12*o105*o113*o29*x+12*o103*o114*o29*x))+6*o153*(x**5)+96*(x**7)))

    
def fprime2(x):
    o45=(r**-2)
    o46=x*x
    o47=o45*o46
    o29=r*r
    o44=m1*m1
    o54=m2*m2
    o103=o29*o54
    o104=o103+o46
    o60=o29*o44
    o61=o46+o60
    o62=(o61**-2)
    o105=(o104**-2)
    o52=o44+o47
    o58=o47+o54
    o123=sqrt(o52)
    o124=sqrt(o58)
    o130=(m1**4)
    o158=(r**4)
    o131=(m2**4)
    o8=(r**-5)
    o35=1/sqrt(o29)
    o114=1/o61
    o108=1/o104
    o53=1/sqrt(o52)
    o59=1/sqrt(o58)
    o122=(o29**1.5)
    o126=o123+o124
    o205=(x**3)
    o147=(x**6)
    o135=(r**6)
    o136=2*o123*o124
    o137=o136+o44+o54
    o152=18*o29*o44
    o153=18*o29*o54
    o154=12*o123*o124*o29
    o155=o152+o153+o154
    o157=(x**4)
    o159=4*o130*o158
    o160=4*o131*o158
    o161=12*o123*o124*o158*o54
    o163=27*o158*o54
    o164=12*o123*o124*o158
    o165=o163+o164
    o166=o165*o44
    o167=o159+o160+o161+o166
    o127=-(o104*o122*o126*o46*o61)
    o129=(r**8)
    o132=o129*o130*o131
    o133=(x**8)
    o134=12*o133
    o139=6*o135*o137*o44*o46*o54
    o156=o147*o155
    o168=o157*o167
    o169=o132+o134+o139+o156+o168
    o170=G*o169
    o171=o127+o170
    o237=(o58**-1.5)
    o240=(o52**-1.5)
    o198=o45*o53*x
    o199=o45*o59*x
    o203=o198+o199
    o204=-(o104*o122*o203*o46*o61)
    o206=-2*o122*o126*o205*o61
    o207=-2*o104*o122*o126*o205
    o208=-2*o104*o122*o126*o61*x
    o210=(x**7)
    o211=96*o210
    o212=12*o123*o59*x
    o213=12*o124*o53*x
    o214=o212+o213
    o215=o147*o214
    o216=2*o123*o45*o59*x
    o217=2*o124*o45*o53*x
    o218=o216+o217
    o219=6*o135*o218*o44*o46*o54
    o220=12*o135*o137*o44*o54*x
    o221=(x**5)
    o222=6*o155*o221
    o223=12*o123*o29*o54*o59*x
    o224=12*o124*o29*o53*o54*x
    o225=12*o123*o29*o59*x
    o226=12*o124*o29*o53*x
    o227=o225+o226
    o228=o227*o44
    o229=o223+o224+o228
    o230=o157*o229
    o231=4*o167*o205
    o232=o211+o215+o219+o220+o222+o230+o231
    o233=G*o232
    o235=o204+o206+o207+o208+o233
    o249=(r**-4)
    o270=-(o240*o249*o46)
    o271=o45*o53
    o272=-(o237*o249*o46)
    o273=o45*o59
    o295=2*o249*o46*o53*o59
    o296=o270+o271
    o297=o124*o296
    o298=o272+o273
    o300=o123*o298
    o301=o295+o297+o300
    
    return 2*(o235*o53*o59-o171*o237*o45*o53*x-o171*o240*o45*o59*x)*(-2*o105*o114*o35*o8*x-2*o108*o35*o62*o8*x)+o108*o114*o35*o8*(2*o235*(-(o237*o45*o53*x)-o240*o45*o59*x)+o53*o59*(-2*o104*o122*o126*o61-4*o122*x*(o104*o203*o61+2*o104*o126*x+2*o126*o61*x)+G*(672*o147+30*o155*o157+12*o214*o221+8*o205*o229+12*o147*o29*o301+12*o167*o46+o157*(12*o158*o301*o44+12*o158*o301*o54)+6*o44*o54*(2*o135*o137+o135*o46*(2*o124*o296+2*o123*o298+4*o249*o46*o53*o59)+4*o135*o218*x))-o122*o46*(o104*(o270+o271+o272+o273)*o61+o126*(2*o104+8*o46+2*o61)+2*o203*(2*o104*x+2*o61*x)))+o171*(2*o237*o240*o249*o46+o59*(-(o240*o45)+3*o249*o46*(o52**-2.5))+o53*(-(o237*o45)+3*o249*o46*(o58**-2.5))))+o171*o35*o53*o59*o8*(8*o105*o46*o62+o114*(-2*o105+8*o46*(o104**-3))+o108*(-2*o62+8*o46*(o61**-3)))
    
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
    
px1 = -pphi/(r)*sin(phi)
py1 = pphi/(r)*cos(phi)
pz1 = 0
px2 = pphi/(r)*sin(phi)
py2 = -pphi/(r)*cos(phi)
pz2 = 0

CX = (m1*x1+m2*x2)/M
CY = (m1*y1+m2*y2)/M
CZ = (m1*z1+m2*z2)/M

CPX = (px1+px2)/M
CPY = (py1+py2)/M
CPZ = (pz1+pz2)/M

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

s0 = str(G) + " " + str(M) + " " + str(L) + " " + str(T) + " 0 0 0\n"
s1 = str(m1) + " " + str(x1) + " " + str(y1) + " " + str(z1) + " " + str(px1) + " " +  str(py1) + " " + str(pz1) + "\n"
s2 = str(m2) + " " + str(x2) + " " + str(y2) + " " + str(z2) + " " + str(px2) + " " +  str(py2) + " " + str(pz2) + "\n"

f.write(s0)
f.write(s1)
f.write(s2)

f.close()

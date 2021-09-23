#----------------------------------------------------------------------------
#               Chaos IC generator 
# File name:    run_nbody.py
# Creator:      Taylor Hugh Morgan, David Neilsen, Jared Jay
# Description:  Generating the Initial Conditions for a choas simulation
#               if LTRACE is set to false then it will automatically run nbodyPN 
#----------------------------------------------------------------------------
import argparse
import sys
from math import *
import subprocess
from scipy import optimize

LTRACE = True

#----------------------------------------------------------------------
#  Set up commandline arguments
#----------------------------------------------------------------------
parser = argparse.ArgumentParser()

parser.add_argument('-t','--time', type=float, default=10.0,metavar='NUM',help="integration time")
parser.add_argument('-d','--time_step', type=float, default=0.1,metavar='NUM',help="time step for integrator")
parser.add_argument('-o','--order', type=str, default="1",
     help='Post-Minkowskian order of the equations',
     choices=["0","1"])
parser.add_argument('--m1', type=float, default=1.0,
       help='Mass of star 1 (Default 1.0)', metavar='NUM')
parser.add_argument('--m2', type=float, default=1.0,
       help='Mass of star 2 (Default 1.0)', metavar='NUM')
parser.add_argument('--m3', type=float, default=1.0,
       help='Mass of star 3 (Default 1.0)', metavar='NUM')
parser.add_argument('-e','--tolerance', type=float, default=1.0e-13,metavar='NUM',help="tolerance for LSODA integrator")
parser.add_argument('-r', '--radius', type=float, default=20.0,
     help='seperation of the binary star system(default 20.0)', metavar = 'NUM')
parser.add_argument('-n','--run_number', type=int, default=0,metavar='INT',help="numerical tag for output files")
parser.add_argument('-f', '--no_com_frame', action="store_true",
     help='Do not use the center of mass frame.')
parser.add_argument('--two_body', action="store_true", help='This a two-body system.')
parser.add_argument('Phi', type=float, default=0.0, help='Phase angle of the binary star')
parser.add_argument('Rho', type=float, default=0.0, help='Impact Parameter')
parser.add_argument('--vc', type=float, default=0.5,metavar='NUM',help="the fraction of the critical velocity that we wish to give to the third object (default=0.5)")

args = parser.parse_args()
r = args.radius
phi = args.Phi
rho = args.Rho
m1 = args.m1
m2 = args.m2
m3 = args.m3


if (args.two_body or args.no_com_frame):
    COM_frame = False
else:
    COM_frame = True

#constants and initial conditions
G = 1.0
c = 1.0
msun = 1.0

m1 = m1*msun
m2 = m2*msun

if args.two_body:
    m3 = 1.0e-15
else:
    m3 = m3*msun

M = m1 + m2
MT = m1 + m2 + m3
mu = m1*m2/(m1+m2)
MU = 2.0*m1*m3/(2.0*m1*m3)
nu = mu/M
pr = 0.0

if args.two_body:
    v3 = 0.0
else:
    v3 = sqrt(G*m1*m2*MT/(2.0*r*m3*(m1+m2)))


#initial position for the object
x1 = r*m2/(m1+m2)*cos(phi)
y1 = r*m2/(m1+m2)*sin(phi)
z1 = 0.0
x2 = -r*m1/(m1+m2)*cos(phi)
y2 = -r*m1/(m1+m2)*sin(phi)
z2 = 0.0

if args.two_body:
    x3 = 1.0e40
else:
    x3 = 60*r

y3 = rho
z3 = 0.0

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

# Initial momentums
px1 = -pphi/(r)*sin(phi)
py1 = pphi/(r)*cos(phi)
pz1 = 0
px2 = pphi/(r)*sin(phi)
py2 = -pphi/(r)*cos(phi)
pz2 = 0
px3 = m3*v3*(-1.0*args.vc)
py3 = 0
pz3 = 0

#Transform to center of momentum frame
if COM_frame == True:
    CX = (m1*x1+m2*x2+m3*x3)/MT
    CY = (m1*y1+m2*y2+m3*y3)/MT
    CZ = (m1*z1+m2*z2+m3*z3)/MT

    CPX = (px1+px2+px3)/MT
    CPY = (py1+py2+py3)/MT
    CPZ = (pz1+pz2+pz3)/MT
else:
    CX = 0.0
    CY = 0.0
    CZ = 0.0

    CPX = 0.0
    CPY = 0.0
    CPZ = 0.0

#position transformation
x1 = x1-CX
y1 = y1-CY
z1 = z1-CZ
x2 = x2-CX
y2 = y2-CY
z2 = z2-CZ
x3 = x3-CX
y3 = y3-CY
z3 = z3-CZ


#momentum transformation
px1 = px1-CPX
py1 = py1-CPY
pz1 = pz1-CPZ
px2 = px2-CPX
py2 = py2-CPY
pz2 = pz2-CPZ
px3 = px3-CPX
py3 = py3-CPY
pz3 = pz3-CPZ

filename = 'ICfile' + str(args.run_number)
#returning the outputs
f = open(filename,'w')
l1 = [m1, x1, y1, z1, px1, py1, pz1]

s1 = str(m1) + " " + str(x1) + " " + str(y1) + " " + str(z1) + " " + str(px1) + " " +  str(py1) + " " + str(pz1) + "\n"
s2 = str(m2) + " " + str(x2) + " " + str(y2) + " " + str(z2) + " " + str(px2) + " " +  str(py2) + " " + str(pz2) + "\n"
s3 = str(m3) + " " + str(x3) + " " + str(y3) + " " + str(z3) + " " + str(px3) + " " +  str(py3) + " " + str(pz3) + "\n"

f.write(s1)
f.write(s2)
f.write(s3)

# Need to close filename so that it is written before it is read by nbodyPN
f.close()

#------------Now ready to run nbodyPN------

if LTRACE:
    print ("Commandline: nbodyPM " + str(args.time) + " " + str(args.time_step) + " " + str(order) + " " + str(args.run_number) + " " + str(args.tolerance) + " " +  str(G) + " " +  str(rho) + " " + str(phi))
else:
    myinput = open(filename)
    myoutput = open('run_nbody.stdout','w')
    p = subprocess.Popen(['nbodyPN',str(args.time),str(args.time_step), 
                      str(order), str(args.run_number),
                      str(args.tolerance), str(G), str(rho), str(phi), str(1)],
                      stdin=myinput, stdout=myoutput,stderr=subprocess.STDOUT)
    p.communicate()  # communicate makes python wait for nbodyPN to run. p.wait()
    myoutput.close()
    myinput.close()

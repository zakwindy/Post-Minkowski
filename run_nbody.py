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

LTRACE = True

#----------------------------------------------------------------------
#  Routines to calculate circular binary initial data
#----------------------------------------------------------------------

#----------------- 0th order-------------------------------------------
def circular_orbit_pn0(r, mu, M):
    q  = r/M
    return sqrt(q)

#----------------- 1st order-------------------------------------------
def circular_orbit_pn1(r, mu, M):
    nu = mu / M
    q  = r/M
    o2=q*q
    pt1=7.071067811865475e-1*sqrt(8.e0*o2-3.9000000000000004e1*q-
        q*sqrt(1.553e3+6.4e1*o2-6.5600000000000005e2*q))
    return pt1


#----------------- 2nd order-------------------------------------------
def circular_orbit_pn2(r, mu, M):
    nu = mu / M
    q  = r/M

    o2=q*q
    o10=q**6
    o4=q**3
    o5=-3.39093e5*o4
    o6=q**4
    o7=-2.6546399999999997e5*o6
    o8=q**5
    o9=1.3593600000000001e5*o8
    o11=-2.3552e4*o10
    o12=2.128054929e9*o10
    o13=q**7
    o14=-2.286574152e9*o13
    o15=q**8
    o16=1.250522344e9*o15
    o17=q**9
    o18=-4.1353344e8*o17
    o19=q**10
    o20=8.7023616e7*o19
    o21=q**11
    o22=-1.0846208e7*o21
    o23=q**12
    o24=6.5536e5*o23
    o25=o12+o14+o16+o18+o20+o22+o24
    o26=sqrt(o25)
    o27=4.409081537009721e1*o26
    o28=o11+o27+o5+o7+o9
    o29=o28**(-3.333333333333333e-1)
    pt2=sqrt(1.7777777777777777e0*o2+
        1.1111111111111112e-1*o28**3.333333333333333e-1-
        1.767e3*o2*o29+6.773333333333333e2*o29*o4-
        9.955555555555556e1*o29*o6+1.666666666666667e0*q)
    return pt2
    
#----------------------------------------------------------------------
#  Set up commandline arguments
#----------------------------------------------------------------------
parser = argparse.ArgumentParser()

parser.add_argument('-t','--time', type=float, default=10.0,metavar='NUM',help="integration time")
parser.add_argument('-d','--time_step', type=float, default=0.1,metavar='NUM',help="time step for integrator")
parser.add_argument('-o','--order', type=str, default="2", 
     help='Post-Newtonian order of the equations', 
     choices=["0","1","2","2.5"])
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

if args.order == "0":
    pphi = circular_orbit_pn0(r, mu, M)
    order = 0
elif args.order == "1":
    pphi = circular_orbit_pn1(r, mu, M)
    order = 1
elif args.order == "2":
    pphi = circular_orbit_pn2(r, mu, M)
    order = 2
elif args.order == "2.5":
    pphi = circular_orbit_pn2(r, mu, M)
    order = 25
else:
    print ("Unknown equation order = " + ags.order)
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

filename = 'ICfile' + str(args.run_number) + '.txt'
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
    print ("Commandline: nbodyPN " + str(args.time) + " " + str(args.time_step) + " " + str(order) + " " + str(args.run_number) + " " + str(args.tolerance) + " " +  str(G) + " " +  str(rho) + " " + str(phi))
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

import argparse
import numpy as np
import math
from scipy import optimize

#----------Set SI units to be used for calculating real values------#
G = 6.647*10**-11
YEAR = 31536000
AU = 1.496 * 10**11
KM = 1000
SOLAR_MASS = 1.989 * 10**30
G_SCALED = 1

M = SOLAR_MASS

#--------Take in command line arguments---#
parser = argparse.ArgumentParser(description = 'Initial parameters for black hole position and mass.')
parser.add_argument('file', type=argparse.FileType('w'), help='name of the file for this program to write to', metavar='FILE')
parser.add_argument('m1', type=float, help='mass of the first black hole', metavar='NUM')
parser.add_argument('m2', type=float, help='mass of the second black hole', metavar='NUM')
parser.add_argument('separation', type=float, help='separation of the binary star system', metavar='SEPARATION')
parser.add_argument('distance_unit', type=str, help='units for the separation of the bodies, examples are AU, KM, M', metavar = 'UNIT')
parser.add_argument('eccentricity', type=float, help='eccentricity of the orbit, must be between 0 and 1', metavar = 'ECC')
parser.add_argument('relativity', type=int, help='whether or not to use PM corrections, 1 for true or 0 for false', metavar = 'REL')

args = parser.parse_args()

mass1 = args.m1
mass2 = args.m2
D = args.separation
unit = args.distance_unit
ecc = args.eccentricity
rela = args.relativity

#----------Set factors for changing parameters based on eccentricity-#
if not (0 <= ecc <= 1) : ecc = 0
qfactor = 1 - ecc
pfactor = math.sqrt((1 + ecc)/(1 - ecc))

#---------Set the units for distance and time---#
if unit == "AU": L = AU
elif unit == "KM": L = KM
elif unit == "M": L = 1

T = math.sqrt(L**3 / (G * M))

#--------Define the functions to solve for initial momentum in PM---#
def f(x):
    o4=(D**-2)
    o5=x*x
    o6=o4*o5
    o3=mass1*mass1
    o10=mass2*mass2
    o14=D*D
    o15=o14*o3
    o16=o15+o5
    o18=o10*o3
    o19=(x**4)
    o20=-2*o19
    o21=o18+o20
    return (-3*G_SCALED*o21*D)/o16+2*G_SCALED*o21*o5*D*(o16**-2)-1/math.sqrt(o10+o6)-1/math.sqrt(o3+o6)

def f1(x):
    o4=mass1*mass1
    o5=D*D
    o6=o4*o5
    o7=x*x
    o8=o6+o7
    o14=(D**-2)
    o15=o14*o7
    o11=(x**3)
    o19=mass2*mass2
    o9=(o8**-2)
    o24=o19*o4
    o25=(x**4)
    o26=-2*o25
    o27=o24+o26
    return (24*G_SCALED*o11*D)/o8+10*G_SCALED*o27*o9*D*x+o14*x*((o15+o19)**-1.5)+o14*x*((o15+o4)**-1.5)-8*G_SCALED*o11*o27*D*(o8**-3)-16*G_SCALED*o9*D*(x**5)
    
def f2(x):
    o4=mass1*mass1
    o5=D*D
    o6=o4*o5
    o7=x*x
    o8=o6+o7
    o17=(D**-2)
    o18=o17*o7
    o19=o18+o4
    o16=(D**-4)
    o24=mass2*mass2
    o25=o18+o24
    o11=(x**4)
    o9=(o8**-3)
    o31=o24*o4
    o32=-2*o11
    o33=o31+o32
    o12=(o8**-2)
    return -208*G_SCALED*o11*o12*D+10*G_SCALED*o12*o33*D+(72*G_SCALED*o7*D)/o8-64*G_SCALED*o33*o7*o9*D-3*o16*o7*(o19**-2.5)+o17*(o19**-1.5)-3*o16*o7*(o25**-2.5)+o17*(o25**-1.5)+48*G_SCALED*o11*o33*D*(o8**-4)+128*G_SCALED*o9*D*(x**6)

#-----------Establish reduced mass and a Newtonian estimate for inital momentum--#
reduced_mass = mass1 * mass2 / ( mass1 + mass2 )
newton_p_theta = ( G_SCALED * mass1 * mass2 * reduced_mass * D )**(0.5)

if rela == 1: p_theta = optimize.newton(f, newton_p_theta, fprime = None, args=(), tol=1e-08, maxiter=50, fprime2 = None, x1=None, rtol=0.0, full_output=False, disp=True)
elif rela == 0: p_theta = newton_p_theta

#-----------Set x and y coordinates and momenta----#
x2, y2 = D / (1 + ( mass2 / mass1 )), 0
x1, y1 = x2 - D, 0
x2 *= qfactor
x1 *= qfactor
px2, py2 = 0, p_theta / D
px1, py1 = 0, -py2
py2 *= pfactor
py1 *= pfactor

line1 = str(G_SCALED) + ' ' + str(M) + ' ' + str(L) + ' ' + str(T) 
line2 = str(mass1) + ' ' + str(x1) + ' ' + str(y1) + ' ' + str(px1) + ' ' + str(py1)
line3 = str(mass2) + ' ' + str(x2) + ' ' + str(y2) + ' ' + str(px2) + ' ' + str(py2)

f = args.file 
f.write(line1)
f.write('\n')
f.write(line2)
f.write('\n')
f.write(line3)
f.close()





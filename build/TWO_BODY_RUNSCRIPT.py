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
    o3=mass1*mass1
    o10=mass2*mass2
    o6=x*x
    o4=D*D
    o5=o3*o4
    o7=o5+o6
    o24=D**-2
    o25=o24*o6
    return G_SCALED*(o7**-2)*(-(o10*o3*o6*(D**6))-3*o10*(mass1**4)*(D**8)+14*o3*(D**4)*(x**4)+10*o4*(x**6))+(D**7)*(-(math.sqrt(o10+o25)/(o10*o4+o6))-math.sqrt(o25+o3)/o7)

def f1(x):
    o3=mass1*mass1
    o4=D*D
    o5=o3*o4
    o6=x*x
    o7=o5+o6
    o9=mass2*mass2
    o10=D**6
    o12=D**4
    o31=D**-2
    o8=o7**-2
    o33=o31*o6
    o34=o3+o33
    o39=o4*o9
    o40=o39+o6
    o42=o33+o9
    return G_SCALED*o8*(-2*o10*o3*o9*x+56*o12*o3*(x**3)+60*o4*(x**5))-4*G_SCALED*x*(o7**-3)*(-(o10*o3*o6*o9)-3*o9*(mass1**4)*(D**8)+14*o12*o3*(x**4)+10*o4*(x**6))+(D**7)*((-1*o31*x)/(o7*math.sqrt(o34))+2*o8*x*math.sqrt(o34)-(1*o31*x)/(o40*math.sqrt(o42))+2*x*(o40**-2)*math.sqrt(o42))
    
def f2(x):
    o3=mass1*mass1
    o6=x*x
    o4=D*D
    o5=o3*o4
    o7=o5+o6
    o9=mass2*mass2
    o10=(D**6)
    o12=(D**4)
    o14=(x**4)
    o18=(o7**-3)
    o27=(mass1**4)
    o28=(D**8)
    o29=-3*o27*o28*o9
    o30=-(o10*o3*o6*o9)
    o31=14*o12*o14*o3
    o32=(x**6)
    o33=10*o32*o4
    o34=o29+o30+o31+o33
    o40=(D**-2)
    o8=(o7**-2)
    o41=o40*o6
    o42=o3+o41
    o39=1/o7
    o45=1/math.sqrt(o42)
    o48=math.sqrt(o42)
    o38=(D**-4)
    o51=o4*o9
    o52=o51+o6
    o54=o41+o9
    o53=1/o52
    o58=1/math.sqrt(o54)
    o57=(o52**-2)
    o62=math.sqrt(o54)
    return -4*G_SCALED*o18*o34+G_SCALED*o8*(300*o14*o4+168*o12*o3*o6-2*o10*o3*o9)+24*G_SCALED*o34*o6*(o7**-4)+(-(o39*o40*o45)-o40*o53*o58-8*o18*o48*o6+4*o40*o57*o58*o6+2*o57*o62+2*o48*o8+4*o40*o45*o6*o8+o38*o39*o6*(o42**-1.5)-8*o6*o62*(o52**-3)+o38*o53*o6*(o54**-1.5))*(D**7)-8*G_SCALED*o18*x*(-2*o10*o3*o9*x+56*o12*o3*(x**3)+60*o4*(x**5))

#-----------Establish reduced mass and a Newtonian estimate for inital momentum--#
reduced_mass = mass1 * mass2 / ( mass1 + mass2 )
newton_p_theta = ( G_SCALED * mass1 * mass2 * reduced_mass * D )**(0.5)

if rela == 1: p_theta = optimize.newton(f, newton_p_theta, fprime = f1, args=(), tol=1e-08, maxiter=50, fprime2 = f2, x1=None, rtol=0.0, full_output=False, disp=True)
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





import argparse
import numpy as np
import math
from scipy import optimize

#----------Set SI units to be used for calculating real values------#
G = 6.647*10**-11
C = 3.00*10**8
YEAR = 31536000
AU = 1.496 * 10**11
KM = 1000
SOLAR_MASS = 1.989 * 10**30
G_SCALED = 1 / (16 * math.pi)

M = SOLAR_MASS

#--------Take in command line arguments---#
parser = argparse.ArgumentParser(description = 'Initial parameters for black hole position and mass.')
parser.add_argument('file', type=argparse.FileType('w'), help='name of the file for this program to write to', metavar='FILE')
parser.add_argument('m1', type=float, help='mass of the first black hole', metavar='NUM')
parser.add_argument('m2', type=float, help='mass of the second black hole', metavar='NUM')
parser.add_argument('separation', type=float, help='separation of the binary star system', metavar='SEPARATION')
parser.add_argument('eccentricity', type=float, help='eccentricity of the orbit, must be between 0 and 1', metavar = 'ECC')
parser.add_argument('relativity', type=int, help='whether or not to use PM corrections, 1 for true or 0 for false', metavar = 'REL')

args = parser.parse_args()

mass1 = args.m1
mass2 = args.m2
D = args.separation
ecc = args.eccentricity
rela = args.relativity

#----------Set factors for changing parameters based on eccentricity-#
if not (0 <= ecc <= 1) : ecc = 0
qfactor = 1 - ecc
pfactor = math.sqrt((1 + ecc)/(1 - ecc))

#---------Set the units for distance and time---#
T = G**2 * M / (C**3)
L = C * T

#--------Define the functions to solve for initial momentum in PM---#

def f(x):
    o29=D*D
    o46=x*x
    o44=mass1*mass1
    o53=mass2*mass2
    o60=(D**-2)
    o61=o46*o60
    o45=o29*o44
    o47=o45+o46
    o54=o29*o53
    o58=o46+o54
    o62=o44+o61
    o104=o53+o61
    o108=math.sqrt(o62)
    o109=math.sqrt(o104)
    o113=(mass1**4)
    o132=(D**4)
    o114=(mass2**4)
    
    return (1*(D**-5)*(-((o108+o109)*o46*o47*o58*(o29**1.5))+G_SCALED*(6*o44*o46*o53*(2*o108*o109+o44+o53)*(D**6)+o113*o114*(D**8)+(4*o113*o132+4*o114*o132+12*o108*o109*o132*o53+o44*(12*o108*o109*o132+27*o132*o53))*(x**4)+(12*o108*o109*o29+18*o29*o44+18*o29*o53)*(x**6)+12*(x**8))))/(o47*o58*math.sqrt(o104)*math.sqrt(o29)*math.sqrt(o62))
    
def fprime(x):
    o29=D*D
    o46=x*x
    o44=mass1*mass1
    o53=mass2*mass2
    o60=(D**-2)
    o61=o46*o60
    o45=o29*o44
    o47=o45+o46
    o54=o29*o53
    o58=o46+o54
    o62=o44+o61
    o103=1/math.sqrt(o62)
    o104=o53+o61
    o105=1/math.sqrt(o104)
    o106=(o29**1.5)
    o112=(x**3)
    o113=math.sqrt(o62)
    o114=math.sqrt(o104)
    o115=o113+o114
    o129=(D**6)
    o165=(D**4)
    o35=1/math.sqrt(o29)
    o52=1/o47
    o59=1/o58
    o166=(mass1**4)
    o168=(mass2**4)
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
    o206=(D**-7)
    o208=-(o106*o115*o46*o47*o58)
    o210=(D**8)
    o211=o166*o168*o210
    o212=(x**8)
    o213=12*o212
    o214=6*o129*o135*o44*o46*o53
    o215=o122*o153
    o216=o155*o197
    o217=o211+o213+o214+o215+o216
    o218=G_SCALED*o217
    o219=o208+o218
    o8=(D**-5)
    
    return -(o103*o206*o219*o35*o52*o59*x*(o104**-1.5))-2*o103*o105*o219*o35*o59*o8*x*(o47**-2)-2*o103*o105*o219*o35*o52*o8*x*(o58**-2)-o105*o206*o219*o35*o52*o59*x*(o62**-1.5)+o103*o105*o35*o52*o59*o8*(-2*o106*o112*o115*o47-2*o106*o112*o115*o58-2*o106*o115*o47*o58*x-o106*o46*o47*o58*(o103*o60*x+o105*o60*x)+G_SCALED*(4*o112*o197+12*o129*o135*o44*o53*x+o122*(12*o105*o113*x+12*o103*o114*x)+6*o129*o44*o46*o53*(2*o105*o113*o60*x+2*o103*o114*o60*x)+o155*(12*o105*o113*o29*o53*x+12*o103*o114*o29*o53*x+o44*(12*o105*o113*o29*x+12*o103*o114*o29*x))+6*o153*(x**5)+96*(x**7)))

    
def fprime2(x):
    o45=(D**-2)
    o46=x*x
    o47=o45*o46
    o29=D*D
    o44=mass1*mass1
    o54=mass2*mass2
    o103=o29*o54
    o104=o103+o46
    o60=o29*o44
    o61=o46+o60
    o62=(o61**-2)
    o105=(o104**-2)
    o52=o44+o47
    o58=o47+o54
    o123=math.sqrt(o52)
    o124=math.sqrt(o58)
    o130=(mass1**4)
    o158=(D**4)
    o131=(mass2**4)
    o8=(D**-5)
    o35=1/math.sqrt(o29)
    o114=1/o61
    o108=1/o104
    o53=1/math.sqrt(o52)
    o59=1/math.sqrt(o58)
    o122=(o29**1.5)
    o126=o123+o124
    o205=(x**3)
    o147=(x**6)
    o135=(D**6)
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
    o129=(D**8)
    o132=o129*o130*o131
    o133=(x**8)
    o134=12*o133
    o139=6*o135*o137*o44*o46*o54
    o156=o147*o155
    o168=o157*o167
    o169=o132+o134+o139+o156+o168
    o170=G_SCALED*o169
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
    o233=G_SCALED*o232
    o235=o204+o206+o207+o208+o233
    o249=(D**-4)
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
    
    return 2*(o235*o53*o59-o171*o237*o45*o53*x-o171*o240*o45*o59*x)*(-2*o105*o114*o35*o8*x-2*o108*o35*o62*o8*x)+o108*o114*o35*o8*(2*o235*(-(o237*o45*o53*x)-o240*o45*o59*x)+o53*o59*(-2*o104*o122*o126*o61-4*o122*x*(o104*o203*o61+2*o104*o126*x+2*o126*o61*x)+G_SCALED*(672*o147+30*o155*o157+12*o214*o221+8*o205*o229+12*o147*o29*o301+12*o167*o46+o157*(12*o158*o301*o44+12*o158*o301*o54)+6*o44*o54*(2*o135*o137+o135*o46*(2*o124*o296+2*o123*o298+4*o249*o46*o53*o59)+4*o135*o218*x))-o122*o46*(o104*(o270+o271+o272+o273)*o61+o126*(2*o104+8*o46+2*o61)+2*o203*(2*o104*x+2*o61*x)))+o171*(2*o237*o240*o249*o46+o59*(-(o240*o45)+3*o249*o46*(o52**-2.5))+o53*(-(o237*o45)+3*o249*o46*(o58**-2.5))))+o171*o35*o53*o59*o8*(8*o105*o46*o62+o114*(-2*o105+8*o46*(o104**-3))+o108*(-2*o62+8*o46*(o61**-3)))

MTOTAL = mass1 + mass1
reduced_mass = mass1 * mass2 / (MTOTAL) 

#-----------Establish reduced mass and a Newtonian estimate for inital momentum--#
newton_p_theta = ( G_SCALED * mass1 * mass2 * reduced_mass * D )**(0.5)

sol = optimize.root_scalar(f, x0=newton_p_theta, fprime=fprime, fprime2=fprime2, method='newton')

pm_p_theta = sol.root

if rela == 1: p_theta = pm_p_theta
elif rela == 0: p_theta = newton_p_theta

#-----------Set x and y coordinates and momenta----#
x2, y2, z2 = D / (1 + ( mass2 / mass1 )), 0, 0
x1, y1, z1 = x2 - D, 0, 0
x2 *= qfactor
x1 *= qfactor
px2, py2, pz2 = 0.0, p_theta / D, 0.0
px1, py1, pz1 = 0.0, -py2, 0.0
py2 *= pfactor
py1 *= pfactor

line1 = str(G_SCALED) + ' ' + str(M) + ' ' + str(L) + ' ' + str(T) 
line2 = str(mass1) + ' ' + str(x1) + ' ' + str(y1) + ' ' + str(z1) + ' '  + str(px1) + ' ' + str(py1) + ' ' + str(pz1)
line3 = str(mass2) + ' ' + str(x2) + ' ' + str(y2) + ' ' + str(z2) + ' ' + str(px2) + ' ' + str(py2) + ' ' + str(pz2)

f = args.file 
f.write(line1)
f.write('\n')
f.write(line2)
f.write('\n')
f.write(line3)
f.close()





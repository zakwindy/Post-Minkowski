#--------------------------------------------------------------------
#               Three body analysis
# File name:    analysis.py
# Creator:      Zackary Windham
# Description:  Analyzes the eccentricity and semi-major axes of a 
#               hierarchical three body system.
#--------------------------------------------------------------------

from math import *
from numpy import *
import argparse
import sys
import pandas
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick

#--------------------------------------------------
#   Read in the file to get data
#--------------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument('File', help='The file that the script reads.')
parser.add_argument('Bodies', default='PNN', help='This describes what the three bodies are. Valid options are PNN (Pulsar and two neutron stars), PNB (a black hole for the third body), PNIB (an intermediate mass black hole for the third body), PNSB (a supermassive black hole for the third body), PBB, PBIB, PBSB, PIBIB, and PIBSB.')
parser.add_argument('Conditions', default='ICL', help='This describes initial conditions. Valid options are ICL (initially circular libration), ICR (rotation), IEL (initially eccentric libration), and IER.')

args = parser.parse_args()
bodies = args.Bodies
conditions = args.Conditions

with open(sys.argv[1], newline='') as f:
    data = pandas.read_csv(f)

#--------------------------------------------------
#   these variables are always constant
#--------------------------------------------------

C = 1.0     # set relativistic c
mSun = 1.0  # set mass to be in solar mass units
    #values in CGS units
C_CGS = 2.998e10
G_CGS = 6.674e-8
mSun_CGS = 1.989e33
AU_CGS = 1.496e13
    #all runs must use G = 1.0
G = 1.0
M = mSun_CGS        #units of mass
L = M*(G_CGS/G)*((C/C_CGS)**2)  #units of length
T = L*C/C_CGS       #units of time
#define the outer eccentricity to always be 0
e_out_initial = 0.0
M_anom_in = 0
M_anom_out = 20*pi/180
i_out_initial = 0

m2 = 1.4    #this mass is always a pulsar

#--------------------------------------------------
#   these variables change for each run
#--------------------------------------------------

if bodies == 'PNN':
    m1 = 1.4
    m3 = 1.4
    a_in_initial = 0.01*AU_CGS/L
    a_out_initial = 0.2*AU_CGS/L
elif bodies == 'PNB':
    m1 = 1.4
    m3 = 30
    a_in_initial = 0.01*AU_CGS/L
    a_out_initial = 0.5*AU_CGS/L
elif bodies == 'PNIB':
    m1 = 1.4
    m3 = 1e3
    a_in_initial = 0.01*AU_CGS/L
    a_out_initial = 2.5*AU_CGS/L
elif bodies == 'PNSB':
    m1 = 1.4
    m3 = 1e6
    a_in_initial = 0.01*AU_CGS/L
    a_out_initial = 10.0*AU_CGS/L
elif bodies == 'PBB':
    m1 = 30
    m3 = 30
    a_in_initial = 0.1*AU_CGS/L
    a_out_initial = 1.0*AU_CGS/L
elif bodies == 'PBIB':
    m1 = 30
    m3 = 1e3
    a_in_initial = 0.1*AU_CGS/L
    a_out_initial = 7.0*AU_CGS/L
elif bodies == 'PBSB':
    m1 = 30
    m3 = 1e6
    a_in_initial = 0.1*AU_CGS/L
    a_out_initial = 40.0*AU_CGS/L
elif bodies == 'PIBIB':
    m1 = 1e3
    m3 = 1e3
    a_in_initial = 0.1*AU_CGS/L
    a_out_initial = 1.2*AU_CGS/L
elif bodies == 'PIBSB':
    m1 = 1e3
    m3 = 1e6
    a_in_initial = 0.1*AU_CGS/L
    a_out_initial = 10.0*AU_CGS/L
else:
    print("Invalid options for Bodies argument. Code did not run.")
    exit()
    
if conditions == 'ICL':
    e_in_initial = 0.01
    i_in_initial = 60*pi/180
    w_in_initial = 60*pi/180
elif conditions == 'ICR':
    e_in_initial = 0.01
    i_in_initial = 60*pi/180
    w_in_initial = 30*pi/180
elif conditions == 'IEL':
    e_in_initial = 0.6
    i_in_initial = 53*pi/180
    w_in_initial = 90*pi/180
elif conditions == 'IER':
    e_in_initial = 0.6
    i_in_initial = 45*pi/180
    w_in_initial = 60*pi/180
else:
    print("Invalid options for Conditions argument. Code did not run.")
    exit()

#--------------------------------------------------
#   derive the period of the orbits
#--------------------------------------------------
M_in = m1+m2
M_out = M_in+m3

T_in = 2*pi*sqrt((a_in_initial**3)/(G*M_in))
T_out = 2*pi*sqrt((a_out_initial**3)/(G*M_out))

#--------------------------------------------------
#   perform the post-process analysis
#--------------------------------------------------
mu_in = m1*m2/M_in
mu_out = M_in*m3/M_out

length = len(data['qx1'])

x1dot = data['px1']/m1 + data['px3']/M_in
y1dot = data['py1']/m1 + data['py3']/M_in
z1dot = data['pz1']/m1 + data['pz3']/M_in
x2dot = data['px2']/m2 + data['px3']/M_in
y2dot = data['py2']/m2 + data['py3']/M_in
z2dot = data['pz2']/m2 + data['pz3']/M_in
x3dot = data['px3']/m3
y3dot = data['py3']/m3
z3dot = data['pz3']/m3

r_in_vec = [data['qx1']-data['qx2'], data['qy1']-data['qy2'], data['qz1']-data['qz2']]
COM = [(data['qx1']*m1+data['qx2']*m2)/M_in, (data['qy1']*m1+data['qy2']*m2)/M_in, (data['qz1']*m1+data['qz2']*m2)/M_in]
r_out_vec = [data['qx3']-COM[0], data['qy3']-COM[1], data['qz3']-COM[2]]

r_in = (r_in_vec[0]**2 + r_in_vec[1]**2 + r_in_vec[2]**2)**0.5
r_out = (r_out_vec[0]**2 + r_out_vec[1]**2 + r_out_vec[2]**2)**0.5

p_in = [m1*x1dot, m1*y1dot, m1*z1dot]
p_out = [m3*x3dot, m3*y3dot, m3*z3dot]

rdot_in = [p_in[0] / mu_in, p_in[1] / mu_in, p_in[2] / mu_in]
rdot_out = [p_out[0] / mu_out, p_out[1] / mu_out, p_out[2] / mu_out]

v_in = (rdot_in[0]**2 + rdot_in[1]**2 + rdot_in[2]**2)**0.5
v_out = (rdot_out[0]**2 + rdot_out[1]**2 + rdot_out[2]**2)**0.5

E_in = 0.5*(v_in**2) - G*M_in/r_in
E_out = 0.5*(v_out**2) - G*M_out/r_out

a_in = -0.5*G*M_in/E_in
a_out = -0.5*G*M_out/E_out

rv_in_vec = [r_in_vec[1]*rdot_in[2] - r_in_vec[2]*rdot_in[1], r_in_vec[2]*rdot_in[0] - r_in_vec[0]*rdot_in[2], r_in_vec[0]*rdot_in[1] - r_in_vec[1]*rdot_in[0]]
rv_out_vec = [r_out_vec[1]*rdot_out[2] - r_out_vec[2]*rdot_out[1], r_out_vec[2]*rdot_out[0] - r_out_vec[0]*rdot_out[2], r_out_vec[0]*rdot_out[1] - r_out_vec[1]*rdot_out[0]]

rv_in = (rv_in_vec[0]**2 + rv_in_vec[1]**2 + rv_in_vec[2]**2)**0.5
rv_out = (rv_out_vec[0]**2 + rv_out_vec[1]**2 + rv_out_vec[2]**2)**0.5

i_in = arccos(rv_in_vec[2]/rv_in)
i_out = arccos(rv_out_vec[2]/rv_out)

e_arg_in = (rv_in**2)/(a_in*G*M_in)
e_arg_out = (rv_out**2)/(a_out*G*M_out)
for i in range(length):
    if e_arg_in[i] > 1.0:
        e_arg_in[i] = 1.0
    if e_arg_out[i] > 1.0:
        e_arg_out[i] = 1.0

e_in = sqrt(1 - e_arg_in)
e_out = sqrt(1 - e_arg_out)

nrv_in_vec = [-rv_in_vec[1], rv_in_vec[0], 0]
nrv_out_vec = [-rv_out_vec[1], rv_out_vec[0], 0]
nrv_in = sqrt(nrv_in_vec[0]**2 + nrv_in_vec[1]**2 + nrv_in_vec[2]**2)
nrv_out = sqrt(nrv_out_vec[0]**2 + nrv_out_vec[1]**2 + nrv_out_vec[2]**2)

arg1_in = nrv_in_vec[0]/nrv_in
arg1_out = nrv_out_vec[0]/nrv_out
arg2_in = (a_in*(1 - e_in**2)-r_in)/(e_in*r_in)
arg2_out = (a_out*(1 - e_out**2)-r_out)/(e_out*r_out)
for i in range(length):
    if abs(arg1_in[i]) > 1.0:
        arg1_in[i] = arg1_in[i]/abs(arg1_in[i])
    if abs(arg1_out[i]) > 1.0:
        arg1_out[i] = arg1_out[i]/abs(arg1_out[i])
    if abs(arg2_in[i]) > 1.0:
        arg2_in[i] = arg2_in[i]/abs(arg2_in[i])
    if abs(arg2_out[i]) > 1.0:
        arg2_out[i] = arg2_out[i]/abs(arg2_out[i])
        
Omega_in = arccos(arg1_in)
Omega_out = arccos(arg1_out)

arg3_in = (r_in_vec[0]*cos(Omega_in) + r_in_vec[1]*sin(Omega_in))/r_in
arg3_out = (r_out_vec[0]*cos(Omega_out) + r_out_vec[1]*sin(Omega_out))/r_out

for i in range(length):
    if abs(arg3_in[i]) > 1.0:
        arg3_in[i] = arg3_in[i]/abs(arg3_in[i])
    if abs(arg3_out[i]) > 1.0:
        arg3_out[i] = arg3_out[i]/abs(arg3_out[i])

f_in = arccos(arg2_in)
f_out = arccos(arg2_out)

theta_in = arccos(arg3_in)
theta_out = arccos(arg3_out)

w_in = theta_in - f_in
w_out = theta_out - f_out

#--------------------------------------------------
#   integrate to get rid of artificial oscillations
#--------------------------------------------------
dt = data['timestep'][1] - data['timestep'][0]
t_end = data['timestep'][length-1]
steps_per_orbit = int(floor(T_out/dt))
num_orbit = int(floor(length / steps_per_orbit))
#FIXME get this to work
'''
for i in range(num_orbit):
    begin = steps_per_orbit*i
    end = steps_per_orbit*(i+1)
    x = data['timestep'][begin:end]
    period = x[end-1] - x[begin]
    y1 = e_in[begin:end]
    sum1 = sum(y1*dt)/period
    e_in[begin:end] = sum1
'''


#--------------------------------------------------
#   transform data to physical units
#--------------------------------------------------
t = data['timestep']*T/(3600*24*365)    #put time data in units of years
a_in_AU = a_in * L / AU_CGS             #put semi-major axes in AU
a_out_AU = a_out * L / AU_CGS

#--------------------------------------------------
#   create plots
#--------------------------------------------------
plt.figure()
plt.plot(t, e_in)
plt.xlabel('Time (years)')
plt.ylabel('Inner Orbit Eccentricity')
plt.savefig('e_in_' + bodies + '_' + conditions + '.png')

plt.figure()
plt.plot(t, e_out)
plt.xlabel('Time (years)')
plt.ylabel('Outer Orbit Eccentricity')
plt.savefig('e_out_' + bodies + '_' + conditions + '.png')

plt.figure()
plt.plot(t, a_in_AU)
plt.xlabel('Time (years)')
plt.ylabel('Inner Orbit Semimajor Axis (AU)')
plt.gca().yaxis.set_major_formatter(mtick.FormatStrFormatter('%.6f'))
plt.savefig('a_in_' + bodies + '_' + conditions + '.png')

plt.figure()
plt.plot(t, a_out_AU)
plt.xlabel('Time (years)')
plt.ylabel('Outer Orbit Semimajor Axis (AU)')
plt.gca().yaxis.set_major_formatter(mtick.FormatStrFormatter('%.6f'))
plt.savefig('a_out_' + bodies + '_' + conditions + '.png')

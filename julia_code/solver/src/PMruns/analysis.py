#--------------------------------------------------------------------
#               Three body analysis
# File name:    analysis.py
# Creator:      Zackary Windham
# Description:  Analyzes the eccentricity and semi-major axes of a 
#               hierarchical three body system.
#--------------------------------------------------------------------

from math import *
from numpy import array
import argparse
import sys
import pandas

#--------------------------------------------------
#   Read in the file to get data
#--------------------------------------------------
with open(sys.argv[1], newline='') as f:
    data = pandas.read_csv(f) 

parser = argparse.ArgumentParser()

#--------------------------------------------------
#   these variables are always constant
#--------------------------------------------------

C = 1.0     # set relativistic c
mSun = 1.0  # set mass to be in solar mass units
    #values in CGS units
C_CGS = 2.998e10
G_CGS = 6.674e-8
mSun_CGS = 1.989e33
AU_CGS = 1.496e14
    #all runs must use G = 1.0
G = 1.0
M = mSun_CGS        #units of mass
L = M*(G_CGS/G)*((C/C_CGS)**2)  #units of length
T = L*C/C_CGS       #units of time
#define the outer eccentricity to always be 0
e_out = 0.0
M_anom_in = 0
M_anom_out = 20*pi/180
i_out = 0

m2 = 1.4    #this is constant every loop

#--------------------------------------------------
#   these variables change for each run
#--------------------------------------------------

m1 = 1.4
m3 = 1.4

a_in = 0.01*AU_CGS/L
a_out = 0.2*AU_CGS/L

e_in = 0.01
i_in = 60*pi/180
w_in = 60*pi/180

#--------------------------------------------------
#   derive the period of the orbits
#--------------------------------------------------
M_in = m1+m2
M_out = M_in+m3

T_in = 2*pi*sqrt((a_in**3)/(G*M_in))
T_out = 2*pi*sqrt((a_out**3)/(G*M_out))

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
COM = [data['qx1']*m1+data['qx2']*m2/M_in, data['qy1']*m1+data['qy2']*m2/M_in, data['qz1']*m1+data['qz2']*m2/M_in]
data['qx3']-COM[1]
r_out_vec = [data['qx3']-COM[0], data['qy3']-COM[2], data['qz3']-COM[2]]

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

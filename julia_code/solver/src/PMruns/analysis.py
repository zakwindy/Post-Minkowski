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
import csv

#--------------------------------------------------
#   Read in the file to get data
#--------------------------------------------------
with open(sys.argv[1], newline='') as f:
    filereader = csv.reader(f, delimiter = ' ', quotechar='|')


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

##these variables change for each run

m1 = 1.4
m2 = 1.4
m3 = 1.4

a_in = 0.01*AU_CGS/L
a_out = 0.2*AU_CGS/L

e_in = 0.01
i_in = 60*pi/180
w_in = 60*pi/180

##derive the period of the orbits
M_in = m1+m2
M_out = M_in+m3

T_in = 2*pi*sqrt((a_in**3)/(G*M_in))
T_out = 2*pi*sqrt((a_out**3)/(G*M_out))

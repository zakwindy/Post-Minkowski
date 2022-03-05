import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

C = 1.0    # set relativistic c
mSun = 1.0  # set mass to be in solar mass units
    #values in CGS units
C_CGS = 2.998e10
G_CGS = 6.674e-8
mSun_CGS = 1.989e33
AU_CGS = 1.496e13
KM_CGS = 1e5
    #all runs must use G = 1.0
G = 1.0
M = mSun_CGS        #units of mass
L = M*(G_CGS/G)*((C/C_CGS)**2)  #units of length
T = L*C/C_CGS       #units of time

#newton files
ni1 = pd.read_csv(r'newtonruns/newton_i_PNN_ICL_1.csv')
ni2 = pd.read_csv(r'newtonruns/newton_i_PNN_ICL_2.csv')
ni3 = pd.read_csv(r'newtonruns/newton_i_PNN_ICL_3.csv')
ni4 = pd.read_csv(r'newtonruns/newton_i_PNN_ICL_4.csv')
ni5 = pd.read_csv(r'newtonruns/newton_i_PNN_ICL_5.csv')
ni6 = pd.read_csv(r'newtonruns/newton_i_PNN_ICL_6.csv')
ni7 = pd.read_csv(r'newtonruns/newton_i_PNN_ICL_7.csv')
ni8 = pd.read_csv(r'newtonruns/newton_i_PNN_ICL_8.csv')
ni9 = pd.read_csv(r'newtonruns/newton_i_PNN_ICL_9.csv')

nein1 = pd.read_csv(r'newtonruns/newton_e_in_PNN_ICL_1.csv')
nein2 = pd.read_csv(r'newtonruns/newton_e_in_PNN_ICL_2.csv')
nein3 = pd.read_csv(r'newtonruns/newton_e_in_PNN_ICL_3.csv')
nein4 = pd.read_csv(r'newtonruns/newton_e_in_PNN_ICL_4.csv')
nein5 = pd.read_csv(r'newtonruns/newton_e_in_PNN_ICL_5.csv')
nein6 = pd.read_csv(r'newtonruns/newton_e_in_PNN_ICL_6.csv')
nein7 = pd.read_csv(r'newtonruns/newton_e_in_PNN_ICL_7.csv')
nein8 = pd.read_csv(r'newtonruns/newton_e_in_PNN_ICL_8.csv')
nein9 = pd.read_csv(r'newtonruns/newton_e_in_PNN_ICL_9.csv')

#PM files
pi1 = pd.read_csv(r'PMruns/PM_i_PNN_ICL_1.csv')
pi2 = pd.read_csv(r'PMruns/PM_i_PNN_ICL_2.csv')
pi3 = pd.read_csv(r'PMruns/PM_i_PNN_ICL_3.csv')
pi4 = pd.read_csv(r'PMruns/PM_i_PNN_ICL_4.csv')
pi5 = pd.read_csv(r'PMruns/PM_i_PNN_ICL_5.csv')
pi6 = pd.read_csv(r'PMruns/PM_i_PNN_ICL_6.csv')
pi7 = pd.read_csv(r'PMruns/PM_i_PNN_ICL_7.csv')
pi8 = pd.read_csv(r'PMruns/PM_i_PNN_ICL_8.csv')
pi9 = pd.read_csv(r'PMruns/PM_i_PNN_ICL_9.csv')

pein1 = pd.read_csv(r'PMruns/PM_e_in_PNN_ICL_1.csv')
pein2 = pd.read_csv(r'PMruns/PM_e_in_PNN_ICL_2.csv')
pein3 = pd.read_csv(r'PMruns/PM_e_in_PNN_ICL_3.csv')
pein4 = pd.read_csv(r'PMruns/PM_e_in_PNN_ICL_4.csv')
pein5 = pd.read_csv(r'PMruns/PM_e_in_PNN_ICL_5.csv')
pein6 = pd.read_csv(r'PMruns/PM_e_in_PNN_ICL_6.csv')
pein7 = pd.read_csv(r'PMruns/PM_e_in_PNN_ICL_7.csv')
pein8 = pd.read_csv(r'PMruns/PM_e_in_PNN_ICL_8.csv')
pein9 = pd.read_csv(r'PMruns/PM_e_in_PNN_ICL_9.csv')

#create the figures
fig, axes = plt.subplots(nrows=9, ncols=2)
fig.suptitle('PNN ICL System')

axes[0,0].set_title('Halved 1 Time')
axes[0,0].plot(ni1['time']*T/3600,ni1['i'])
axes[0,0].plot(pi1['time']*T/3600,pi1['i'])
axes[0,0].set_xlabel('Time (hours)')
axes[0,0].set_ylabel('i (degrees)')

axes[0,1].set_title('Halved 1 Time')
axes[0,1].plot(nein1['time']*T/3600,nein1['e_in'])
axes[0,1].plot(pein1['time']*T/3600,pein1['e_in'])
axes[0,1].set_xlabel('Time (hours)')
axes[0,1].set_ylabel('Eccentricity')

axes[1,0].set_title('Halved 2 Times')
axes[1,0].plot(ni2['time']*T/3600,ni2['i'])
axes[1,0].plot(pi2['time']*T/3600,pi2['i'])
axes[1,0].set_xlabel('Time (hours)')
axes[1,0].set_ylabel('i (degrees)')

axes[1,1].set_title('Halved 2 Times')
axes[1,1].plot(nein2['time']*T/3600,nein2['e_in'])
axes[1,1].plot(pein2['time']*T/3600,pein2['e_in'])
axes[1,1].set_xlabel('Time (hours)')
axes[1,1].set_ylabel('Eccentricity')

axes[2,0].set_title('Halved 3 Times')
axes[2,0].plot(ni3['time']*T/3600,ni3['i'])
axes[2,0].plot(pi3['time']*T/3600,pi3['i'])
axes[2,0].set_xlabel('Time (hours)')
axes[2,0].set_ylabel('i (degrees)')

axes[2,1].set_title('Halved 3 Times')
axes[2,1].plot(nein3['time']*T/3600,nein3['e_in'])
axes[2,1].plot(pein3['time']*T/3600,pein3['e_in'])
axes[2,1].set_xlabel('Time (hours)')
axes[2,1].set_ylabel('Eccentricity')

axes[3,0].set_title('Halved 4 Times')
axes[3,0].plot(ni4['time']*T/3600,ni4['i'])
axes[3,0].plot(pi4['time']*T/3600,pi4['i'])
axes[3,0].set_xlabel('Time (hours)')
axes[3,0].set_ylabel('i (degrees)')

axes[3,1].set_title('Halved 4 Times')
axes[3,1].plot(nein4['time']*T/3600,nein4['e_in'])
axes[3,1].plot(pein4['time']*T/3600,pein4['e_in'])
axes[3,1].set_xlabel('Time (hours)')
axes[3,1].set_ylabel('Eccentricity')

axes[4,0].set_title('Halved 5 Times')
axes[4,0].plot(ni5['time']*T/3600,ni5['i'])
axes[4,0].plot(pi5['time']*T/3600,pi5['i'])
axes[4,0].set_xlabel('Time (hours)')
axes[4,0].set_ylabel('i (degrees)')

axes[4,1].set_title('Halved 5 Times')
axes[4,1].plot(nein5['time']*T/3600,nein5['e_in'])
axes[4,1].plot(pein5['time']*T/3600,pein5['e_in'])
axes[4,1].set_xlabel('Time (hours)')
axes[4,1].set_ylabel('Eccentricity')

axes[5,0].set_title('Halved 6 Times')
axes[5,0].plot(ni6['time']*T/3600,ni6['i'])
axes[5,0].plot(pi6['time']*T/3600,pi6['i'])
axes[5,0].set_xlabel('Time (hours)')
axes[5,0].set_ylabel('i (degrees)')

axes[5,1].set_title('Halved 6 Times')
axes[5,1].plot(nein6['time']*T/3600,nein6['e_in'])
axes[5,1].plot(pein6['time']*T/3600,pein6['e_in'])
axes[5,1].set_xlabel('Time (hours)')
axes[5,1].set_ylabel('Eccentricity')

axes[6,0].set_title('Halved 7 Times')
axes[6,0].plot(ni7['time']*T/3600,ni7['i'])
axes[6,0].plot(pi7['time']*T/3600,pi7['i'])
axes[6,0].set_xlabel('Time (hours)')
axes[6,0].set_ylabel('i (degrees)')

axes[6,1].set_title('Halved 7 Times')
axes[6,1].plot(nein7['time']*T/3600,nein7['e_in'])
axes[6,1].plot(pein7['time']*T/3600,pein7['e_in'])
axes[6,1].set_xlabel('Time (hours)')
axes[6,1].set_ylabel('Eccentricity')

axes[7,0].set_title('Halved 8 Times')
axes[7,0].plot(ni8['time']*T/3600,ni8['i'])
axes[7,0].plot(pi8['time']*T/3600,pi8['i'])
axes[7,0].set_xlabel('Time (hours)')
axes[7,0].set_ylabel('i (degrees)')

axes[7,1].set_title('Halved 8 Times')
axes[7,1].plot(nein8['time']*T/3600,nein8['e_in'])
axes[7,1].plot(pein8['time']*T/3600,pein8['e_in'])
axes[7,1].set_xlabel('Time (hours)')
axes[7,1].set_ylabel('Eccentricity')

axes[8,0].set_title('Halved 9 Times')
axes[8,0].plot(ni9['time']*T/3600,ni9['i'])
axes[8,0].plot(pi9['time']*T/3600,pi9['i'])
axes[8,0].set_xlabel('Time (hours)')
axes[8,0].set_ylabel('i (degrees)')

axes[8,1].set_title('Halved 9 Times')
axes[8,1].plot(nein9['time']*T/3600,nein9['e_in'])
axes[8,1].plot(pein9['time']*T/3600,pein9['e_in'])
axes[8,1].set_xlabel('Time (hours)')
axes[8,1].set_ylabel('Eccentricity')

plt.tight_layout()
fig.savefig('PNN_ICL.pdf')

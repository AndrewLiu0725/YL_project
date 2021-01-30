# ===============================================================================
# Copyright 2021 An-Jun Liu
# Last Modified Date: 01/29/2021
# ===============================================================================
import numpy as np 
import matplotlib.pyplot as plt
import sys
import time
import datetime
from RBC_Utilities import calcDoubletFraction, getInstrinsicViscosity, getStress

"""
This cose is to plot phase diagram for two-cell system
"""

# Make mean doublet fraction (ensemble) versus Ca plot here

start_time = time.time()

ncycle = 2000
#phis = [3.8, 4.7, 6.0, 6.9]
vol = 746.3163
phis  =[]
for x in range(30, 41):
    phi = 2*vol/(24*x**2)
    phis.append(float(str(phi*100)[:3]))
Ca_list = [(i+1)*0.01 for i in range(20)]
angles = [90-10*i for i in range(18)]
r = 1

dfa_g = np.zeros((len(phis), len(Ca_list)))

for phi_index, phi in enumerate(phis):
    print('phi = ', phi)
    for Ca_index, Ca in enumerate(Ca_list):
        for angle_index, angle in enumerate(angles):
            result = calcDoubletFraction(phi, Ca, 1, r, 0, angle, 0)
            dfa_g[phi_index, Ca_index] += np.mean(result[0][0][-int(result[1]/2):])
dfa_g = dfa_g/len(angles) # average over the ensemble

plt.figure(figsize = (16,12))
for phi_index, phi in enumerate(phis):
    plt.plot(Ca_list, dfa_g[phi_index, :] ,label = 'phi = {}'.format(phi))
    
plt.xlabel("Ca", fontsize = 30)
plt.xticks(fontsize = 20)
plt.ylabel("Doublet Fraction", fontsize = 30)
plt.yticks(fontsize = 20)
plt.legend(prop={'size': 20})
plt.title("Doublet Fraction vs Ca (Two-cell system)\ncriteria_r = {}Dm, criteria_T = 1t_rot, second half time series".format(r), fontsize = 30)
plt.tight_layout()
plt.savefig("./Pictures/TwoCellSystem_DoubletFraction_vs_Ca_SecondHalf_r_{}.png".format(r), dpi = 300)

from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from matplotlib import rcParams


def Z(x, y, df):
    return df[x, y]

df = dfa_g
text1 = 'Doublet Fraction'
text2 = 'DoubletFraction'

        
fig = plt.figure(figsize = (16,12))
ax = fig.gca(projection='3d')
index_X, index_Y = np.meshgrid(np.array(list(range(len(phis)))), np.array(list(range(len(Ca_list)))))
X, Y = np.meshgrid(np.array(phis), np.array(Ca_list))
surf = ax.plot_surface(Y, X, Z(index_X, index_Y, df), cmap = cm.coolwarm, edgecolor = 'black', linewidth = 2)

ax.set_zlim(0.0, np.max(df))
ax.set_xlabel('Ca', fontsize = 30)
ax.set_ylabel('phi', fontsize = 30)
ax.set_zlabel('Doublet fraction', fontsize = 30)
rcParams['axes.labelpad'] = 20
ax.tick_params(labelsize = 20)
ax.set_title('Phase Diagram of phi vs Ca\n'+text1+'\ncriteria_r = {}Dm, criteria_T = 1t_rot, second half time series'.format(r), fontsize = 30)
ax.view_init(azim = -80, elev = 45)
fig.colorbar(surf, ax = ax)
fig.tight_layout()
plt.savefig('./Pictures/TwoCellSystem_PhaseDiagram_3D_{}_SecondHalf_r_{}.png'.format(text2, r), dpi = 300)


fig, ax = plt.subplots(figsize = (16,12))
grid_phi = np.array(phis+[phis[-1]+1])-0.5
grid_Ca = np.array(Ca_list+[Ca_list[-1]+0.01])-0.005
wider_phis, grid_phi  =[], []
for x in range(29, 42):
    phi = 2*vol/(24*x**2)
    wider_phis.append(float(str(phi*100)[:3]))
for i in range(len(phis)+1):
    grid_phi.append((wider_phis[i]+wider_phis[i+1])/2)
wider_Ca_list = [(i)*0.01 for i in range(22)]
grid_Ca = []
for i in range(len(Ca_list)+1):
    grid_Ca.append((wider_Ca_list[i]+wider_Ca_list[i+1])/2)
grid_phi = np.array(grid_phi)
grid_Ca = np.array(grid_Ca)
		
X, Y = np.meshgrid(grid_phi, grid_Ca)
im = ax.pcolormesh(Y, X, Z(index_X, index_Y, df), cmap = cm.coolwarm)
cb = fig.colorbar(im, ax = ax)
cb.ax.tick_params(labelsize = 25)
ax.set_title('Phase Diagram of phi vs Ca\n'+text1+'\ncriteria_r = {}Dm, criteria_T = 1t_rot, second half time series'.format(r), fontsize = 30)
ax.set_xlabel('Ca', fontsize = 30)
ax.set_ylabel('phi', fontsize = 30)
ax.tick_params(labelsize = 25)
fig.tight_layout()
plt.savefig('./Pictures/TwoCellSystem_PhaseDiagram_2D_{}_SecondHalf_r_{}.png'.format(text2, r), dpi = 300)

print('Total time elapsed = {}'.format(str(datetime.timedelta(seconds=time.time()-start_time))))
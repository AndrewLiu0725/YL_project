# ===============================================================================
# Copyright 2021 An-Jun Liu
# Last Modified Date: 02/04/2021
# ===============================================================================
import numpy as np
import matplotlib.pyplot as plt
import time
import datetime
from RBC_Utilities import calcDoubletFraction
import os
import sys 
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from matplotlib import rcParams



"""
This code is to make two plots: doublet fraction vs Ca plot and its phasediagram (in colormap).
"""

def Z(x, y, df):
    return df[x, y]

start_time = time.time()

path = "/raid6/ctliao/Data/HI_ordering/"

# parser
def parser(string):
    # output [phi, Ca, ensemble id]
    string = string.split("-")
    ensemble_id = int(string[1])
    # format: h24_phi4.4989_Re0.1_Ca0.06_WCA1_zero0.8-8
    if "_" in string[0]:
        string = string[0].split("_")
        phi = float(string[1][3:])
        Ca = float(string[3][2:])
    # format: h24phi4.9488Re0.1Ca0.06WCA1zero0.8-4
    else:
        phi = float(string[0][string[0].find('i')+1:string[0].find('R')])
        Ca = float(string[0][string[0].find('a')+1:string[0].find('W')])
    return [phi, Ca, ensemble_id]


parameter_set = {}
# create parameter_set (two layer dict)
for fn in os.listdir(path):
    if (fn[0] == "h") and (os.path.isdir(path+fn)):
        result = parser(fn)
        [phi, Ca, ensemble_id] = result
        if phi in parameter_set.keys():
            if Ca in parameter_set[phi].keys():
                parameter_set[phi][Ca].append(ensemble_id)
            else:
                parameter_set[phi][Ca] = [ensemble_id]
        else:
            parameter_set[phi] = {}
            parameter_set[phi][Ca] = [ensemble_id]


phis = []
Cas = set([i*0.01 for i in range(1, 21)])
for phi in parameter_set.keys():
    if (len(parameter_set[phi].keys()) > 5):
        phis.append(phi)
        Cas = Cas.intersection(set(list(parameter_set[phi].keys())))
phis.sort()
Cas = list(Cas)
Cas.sort()


# make plots
r = 1
avg_df = np.zeros((len(phis), len(Cas)))
for phi_index, phi in enumerate(phis):
    for Ca_index, Ca in enumerate(Cas):
        ensemble_count = 0
        sum_df = 0
        # run over ensemble
        for ensemble_id in parameter_set[phi][Ca]:
            try:
                result = calcDoubletFraction(phi, Ca, 1, r, 0, ensemble_id, 1)
                sum_df += np.mean(result[0][0][-int(result[1]/2):]) # last half
                ensemble_count += 1

            except KeyboardInterrupt:
                print("Pressed ctrl+C\n")
                sys.exit()

            except:
                pass
                #print("Error: no preprocessed data in (phi = {}, Ca = {}, ensemble_id = {})".format(phi, Ca, ensemble_id))

        if (ensemble_count > 0):
            avg_df[phi_index, Ca_index] = sum_df/ensemble_count
        else:
            print("Error: no data in (phi = {}, Ca = {})".format(phi, Ca))



# doublet fraction vs Ca plot
plt.figure(figsize = (16,12))
for phi_index, phi in enumerate(phis):
    plt.plot(Cas, avg_df[phi_index, :] ,label = 'phi = {}'.format(phi))
    
plt.xlabel("Ca", fontsize = 30)
plt.xticks(fontsize = 20)
plt.ylabel("Doublet Fraction", fontsize = 30)
plt.yticks(fontsize = 20)
plt.legend(prop={'size': 20})
plt.title("Doublet Fraction vs Ca (Suspension system)\ncriteria_r = {}Dm, criteria_T = 1t_rot, second half time series".format(r), fontsize = 30)
plt.tight_layout()
plt.savefig("./Pictures/SuspensionSystem_DoubletFraction_vs_Ca_EnsembleAveraged_SecondHalf_r_{}.png".format(r), dpi = 300)
plt.close()



# 3D phase diagram  
df = avg_df
text1 = 'Doublet Fraction'
text2 = 'DoubletFraction'
      
fig = plt.figure(figsize = (16,12))
ax = fig.gca(projection='3d')
index_X, index_Y = np.meshgrid(np.array(list(range(len(phis)))), np.array(list(range(len(Cas)))))
X, Y = np.meshgrid(np.array(phis), np.array(Cas))
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
plt.savefig('./Pictures/SuspensionSystem_PhaseDiagram_3D_{}_EnsembleAveraged_SecondHalf_r_{}.png'.format(text2, r), dpi = 300)
plt.close()



# 2D phase diagram
fig, ax = plt.subplots(figsize = (16,12))
grid_phi = np.array([phis[0]-(phis[1]-phis[0])/2] + [(phis[i]+phis[i+1])/2 for i in range(len(phis)-1)] + [phis[-1] + (phis[-1]-phis[-2])/2])
grid_Ca = np.array([Cas[0]-(Cas[1]-Cas[0])/2] + [(Cas[i]+Cas[i+1])/2 for i in range(len(Cas)-1)] + [Cas[-1] + (Cas[-1]-Cas[-2])/2])
		
X, Y = np.meshgrid(grid_phi, grid_Ca)
im = ax.pcolormesh(Y, X, Z(index_X, index_Y, df), cmap = cm.coolwarm)
cb = fig.colorbar(im, ax = ax)
cb.ax.tick_params(labelsize = 25)
ax.set_title('Phase Diagram of phi vs Ca\n'+text1+'\ncriteria_r = {}Dm, criteria_T = 1t_rot, second half time series'.format(r), fontsize = 30)
ax.set_xlabel('Ca', fontsize = 30)
ax.set_ylabel('phi', fontsize = 30)
ax.tick_params(labelsize = 25)
fig.tight_layout()
plt.savefig('./Pictures/SuspensionSystem_PhaseDiagram_2D_{}_EnsembleAveraged_SecondHalf_r_{}.png'.format(text2, r), dpi = 300)
plt.close()


print('\nTotal time elapsed = {}'.format(str(datetime.timedelta(seconds=time.time()-start_time))))
# ===============================================================================
# Copyright 2021 An-Jun Liu
# Last Modified Date: 08/16/2022
# ===============================================================================
import numpy as np 
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

"""
This code is to make the plots of intrinsic and relative viscosity vs Ca and phi for two-cell system.
The data is pre-calculated. 
"""

MAKE_PLOT = 1
plt.rcParams['font.family'] = 'DeJavu Serif'
plt.rcParams['font.serif'] = ['Times New Roman']

data = np.load("Data/viscosity_data_0.5.npy")

# Parameters
r = 1.0
vol = 746.3163
phis  =[]
for x in range(30, 41):
    phi = 2*vol/(24*x**2)
    #phis.append(float(str(phi*100)[:3]))
    phis.append(round(phi*100, 1))
phi_range = np.array(phis)
Ca_range = np.array([i*0.01 for i in range(1, 21)])

### vs Ca
for i in range(2): # run over iv and rv
    for j in range(2): # with or wihtout errorbar
        fig, ax = plt.subplots(figsize = (12, 9))
        for phi_index, phi in reversed(list(enumerate(phi_range))):
            if phi_index % 2 == 0:
                continue
            if j == 0:
                ax.errorbar(Ca_range, data[i, 0, phi_index, :], yerr = data[i, 1, phi_index, :], marker = 's', label = r'$\phi$'+' = {}%'.format(phi), capsize = 2)
            else:
                ax.plot(Ca_range, data[i, 0, phi_index, :], marker = 's', label = r'$\phi$'+' = {}%'.format(phi))


        #plt.title("{} vs Ca (two-cell system)".format(r'$\eta _{int}$' if i == 0 else r'$\eta _{rel}$'), fontsize = 30)
        ax.set_xlabel("Ca", fontsize = 30)
        ax.set_ylabel(r'$\left[ \eta \right]$' if i == 0 else r'$\eta _{rel}$', fontsize = 30)
        ax.tick_params(labelsize=22.5)
        ax.tick_params(which='both', width=3)
        ax.xaxis.set_major_locator(MaxNLocator(5))
        ax.xaxis.set_minor_locator(MaxNLocator(10))
        ax.yaxis.set_major_locator(MaxNLocator(5)) 
        ax.yaxis.set_minor_locator(MaxNLocator(10))
        ax.legend(fontsize = 20, frameon = False)
        if MAKE_PLOT:
            fig.savefig("./Pictures/Viscosity_vs_phi_Ca/TwoCellSystem_{}_vs_Ca_{}ErrorBar.png".format("IntrinsicViscosity" if i == 0 else "RelativeViscosity", "with" if j == 0 else "without"), dpi = 300)
            plt.close()
        else:
            plt.show()

### vs phi
for i in range(2): # run over iv and rv
    for j in range(2): # with or wihtout errorbar
        fig, ax = plt.subplots(figsize = (12, 9))
        for Ca_index, Ca in enumerate(Ca_range):
            if not Ca in [0.03, 0.06, 0.08, 0.1, 0.18]: continue
            #if Ca_index % 4 != 0: continue
            if j == 0:
                ax.errorbar(phi_range, data[i, 0, :, Ca_index], yerr = data[i, 1, :, Ca_index], marker = 's', label = 'Ca = {}'.format(Ca), capsize = 2)
            else:
                ax.plot(phi_range, data[i, 0, :, Ca_index], marker = 's', label = 'Ca = {}'.format(Ca))


        #plt.title("{} vs {} (two-cell system)".format(r'$\eta _{int}$' if i == 0 else r'$\eta _{rel}$', r'$\phi$'), fontsize = 30)
        ax.set_xlabel("{} (%)".format(r'$\phi$'), fontsize = 30)
        ax.set_ylabel(r'$\left[ \eta \right]$' if i == 0 else r'$\eta _{rel}$', fontsize = 30)
        ax.tick_params(labelsize=22.5)
        ax.tick_params(which='both', width=3)
        ax.xaxis.set_major_locator(MaxNLocator(5))
        ax.xaxis.set_minor_locator(MaxNLocator(10))
        ax.yaxis.set_major_locator(MaxNLocator(5)) 
        ax.yaxis.set_minor_locator(MaxNLocator(10))
        ax.legend(fontsize = 20, frameon = False)
        if MAKE_PLOT:
            fig.savefig("./Pictures/Viscosity_vs_phi_Ca/TwoCellSystem_{}_vs_phi_{}ErrorBar.png".format("IntrinsicViscosity" if i == 0 else "RelativeViscosity", "with" if j == 0 else "without"), dpi = 300)
            plt.close()
        else:
            plt.show()
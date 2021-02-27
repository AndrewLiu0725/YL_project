# ===============================================================================
# Copyright 2021 An-Jun Liu
# Last Modified Date: 02/21/2021
# ===============================================================================
import numpy as np 
import matplotlib.pyplot as plt

"""
This code is to make the plots of intrinsic and relative viscosity vs Ca and phi for suspension system.
The data is pre-calculated. 
"""

MAKE_PLOT = 1

data = {}
phi_list = [2.3994, 2.9993, 3.4492, 3.8991, 4.9488, 5.9986]

for phi in phi_list:
    data[phi] = np.loadtxt("./Data/data/phi{}.txt".format(str(phi)))

# vs Ca plot
for i in range(2): # iv or rv
    for j in range(2): # with or without errorbar
        plt.figure(figsize=(12,9))

        for phi in phi_list:
            #if phi == 3.4492: continue
            if j == 0:
                plt.errorbar(data[phi][:, 0], data[phi][:, 1 if i == 0 else 3], yerr = data[phi][:, 2 if i == 0 else 4], label = r'$\phi$'+" = {}%".format(phi), capsize = 2)
            else:
                plt.plot(data[phi][:, 0], data[phi][:, 1 if i == 0 else 3], label = r'$\phi$'+" = {}%".format(phi))
        plt.title("{} vs Ca (suspension system)".format(r'$\eta _{int}$' if i == 0 else r'$\eta _{rel}$'), fontsize = 30)
        plt.xlabel("Ca", fontsize = 20)
        plt.xticks(fontsize = 15)
        plt.ylabel(r'$\eta _{int}$' if i == 0 else r'$\eta _{rel}$', fontsize = 20)
        plt.yticks(fontsize = 15)
        plt.legend(fontsize = 15)
        if MAKE_PLOT:
            plt.savefig("./Pictures/Suspension/SuspensionSystem_{}_vs_Ca_{}ErrorBar.png".format("IntrinsicViscosity" if i == 0 else "RelativeViscosity", "with" if j == 0 else "without"), dpi = 300)
            plt.close()
        else:
            plt.show()




Ca_list = data[phi_list[0]][:, 0]
data_Ca = {}
for Ca in Ca_list:
    data_Ca[Ca] = np.zeros((len(phi_list), 5))

for phi_index, phi in enumerate(phi_list):
    for Ca_index, Ca in enumerate(Ca_list):
        data_Ca[Ca][phi_index, 0] = phi
        data_Ca[Ca][phi_index, 1:] = data[phi][Ca_index, 1:]


#  vs phi plots
for i in range(2): # iv or rv
    for j in range(2): # with or without errorbar
        plt.figure(figsize=(12,9))

        for Ca in Ca_list:
            if not Ca in [0.03, 0.08, 0.1, 0.14, 0.18]: continue
            if j == 0:
                plt.errorbar(data_Ca[Ca][:, 0], data_Ca[Ca][:, 1 if i == 0 else 3], yerr = data_Ca[Ca][:, 2 if i == 0 else 4], label = "Ca = {}".format(Ca), capsize = 2)
            else:
                plt.plot(data_Ca[Ca][:, 0], data_Ca[Ca][:, 1 if i == 0 else 3], label = "Ca = {}".format(Ca))
        plt.title("{} vs {} (suspension system)".format(r'$\eta _{int}$' if i == 0 else r'$\eta _{rel}$', r'$\phi$'), fontsize = 30)
        plt.xlabel(r'$\phi$'+"(%)", fontsize = 20)
        plt.xticks(fontsize = 15)
        plt.ylabel(r'$\eta _{int}$' if i == 0 else r'$\eta _{rel}$', fontsize = 20)
        plt.yticks(fontsize = 15)
        plt.legend(fontsize = 15)
        if MAKE_PLOT:
            plt.savefig("./Pictures/Suspension/SuspensionSystem_{}_vs_phi_{}ErrorBar.png".format("IntrinsicViscosity" if i == 0 else "RelativeViscosity", "with" if j == 0 else "without"), dpi = 300)
            plt.close()
        else:
            plt.show()
# ===============================================================================
# Copyright 2021 An-Jun Liu
# Last Modified Date: 08/22/2022
# ===============================================================================
import matplotlib.pyplot as plt 
from matplotlib.ticker import MaxNLocator
import numpy as np 
import pickle

plt.rcParams['font.family'] = 'DeJavu Serif'
plt.rcParams['font.serif'] = ['Times New Roman']

# load data
columns, rows = 2, 3

data = [[], []]
system_list = ['Suspension', 'TwoCell']
variable_list = ['AverageDF', 'IV', 'RV']
for i in range(columns):
    for j in range(rows):
        with open("Data/{}_{}_EnsembleAveraged.pickle".format(variable_list[j], system_list[i]), 'rb') as handle:
            data[i].append(pickle.load(handle))

fig, axs = plt.subplots(rows, columns, figsize = (7.5*columns, 6*rows))

system_list = ["Suspension", "Two-Cell"]
variable_list = ["Doublet Fraction", "Intrinsic Viscosity", "Relative Viscosity"]
symbol_list = [r'$\Phi$', r'$\left[ \eta \right]$', r'$\eta _{rel}$']

slicing = 10

for i in range(columns): # system
    factor = 1 if i else 4000/3669

    if i:
        phis = [5.7, 4.0, 4.0]
        Cas = [0.1, 0.02, 0.17]
    else:
        phis = [2.9993, 2.9993, 5.9986]
        Cas = [0.08, 0.16, 0.08]

    for j in range(rows): # variable
        for k in range(3):
            phi, Ca = phis[k], Cas[k]
            ts = data[i][j][phi][Ca][0]
            axs[j, i].plot((np.arange(len(ts))*factor)[::slicing], ts[::slicing], label = "{} = {}%, Ca = {}".format(r'$\phi$', round(phi,1), Ca))
        axs[j, i].set_xlabel("{}".format(r'$\gamma _{rot}$'), fontsize = 20)
        axs[j, i].set_ylabel(symbol_list[j], fontsize = 20)
        axs[j, i].set_title("{} Time Series ({} System)".format(variable_list[j], system_list[i]), fontsize = 20)
        axs[j, i].tick_params(which='both', labelsize=15, width=2, length=8, direction='in')
        axs[j, i].xaxis.set_major_locator(MaxNLocator(5))
        axs[j, i].xaxis.set_minor_locator(MaxNLocator(10))
        axs[j, i].yaxis.set_major_locator(MaxNLocator(5)) 
        #axs[j, i].yaxis.set_minor_locator(MaxNLocator(10))
        axs[j, i].legend(frameon=False, fontsize=15)
fig.tight_layout()
#plt.show()
fig.savefig("Pictures/Manuscript/fig_S5.png", dpi = 200)
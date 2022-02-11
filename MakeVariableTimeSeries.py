# ===============================================================================
# Copyright 2021 An-Jun Liu
# Last Modified Date: 02/10/2022
# ===============================================================================
import matplotlib.pyplot as plt 
import numpy as np 
import pickle

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
            axs[j, i].plot((np.arange(len(ts))*factor)[::slicing], ts[::slicing], label = "{} = {}, Ca = {}".format(r'$\phi$', round(phi,1), Ca))
        axs[j, i].set_xlabel("{}".format(r'$\dot \gamma t$'), fontsize = 15)
        axs[j, i].set_ylabel(symbol_list[j], fontsize = 15)
        axs[j, i].set_title("{} Time Series ({} System)".format(variable_list[j], system_list[i]), fontsize = 15)
        axs[j, i].tick_params(labelsize = 12)
        axs[j, i].legend()
fig.tight_layout()
#plt.show()
plt.savefig("Pictures/Manuscript/Variables_TimeSeries.png", dpi = 200)
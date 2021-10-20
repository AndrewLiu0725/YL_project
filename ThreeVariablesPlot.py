# ===============================================================================
# Copyright 2021 An-Jun Liu
# Last Modified Date: 10/19/2021
# ===============================================================================
import matplotlib.pyplot as plt 
import numpy as np 
from RBC_Utilities import calcDoubletFraction, getIntrinsicViscosity, getRelativeViscosity
import time
import datetime
import random

start_time = time.time()

### Test part
# list df(phi, Ca, ensemble_id)
r = 1
#parameter_set = [[6.4, 0.02], [6.4, 0.07], [6.4, 0.15], [3.8, 0.07]]
parameter_set = [[6.4, 0.02], [6.4, 0.07], [6.4, 0.15]]

'''
print("===============================================================================")
print("Start listing df(phi, Ca, ensemble_id)")
for [phi, Ca] in parameter_set:
    for angle in [90-10*i for i in range(18)]:
        result = calcDoubletFraction(phi, Ca, 1, r, 0, angle, 0)
        #print("phi = {}, Ca = {}, angle = {}: df = {}".format(phi, Ca, angle, np.mean(result[0][0])))
    #print('\n')
print("Finish listing df(phi, Ca, ensemble_id)")
print("===============================================================================")
'''

# make doublet dominant and singlet dominant plot for each set of phi and Ca (also error bar version)
k1, k2 = 0.2, 0.8 # max avg_df for singlet dominance, min avg_df for doublet dominance
columns, rows = 2, 3
slicing = 40 # for plotting only
st = 1 # to avoid the extreme large initial value for viscosity data
symbol_list = ['df', r'$\eta _{rel}$', r'$\left[ \eta \right]$']

print("===============================================================================")
print("Start finding singlet and doublet dominance")

angleSelected = [[] for _ in range(len(parameter_set))]
timeTrajectories = [[] for _ in range(len(parameter_set))]

for set_count, [phi, Ca] in enumerate(parameter_set):
    fig1, ax1 = plt.subplots(rows, columns, figsize = (7.5*columns, 6*rows)) # for doublet dominance
    fig2, ax2 = plt.subplots(rows, columns, figsize = (7.5*columns, 6*rows)) # for singlet dominance
    numberDoubletDominance, numberSingletDominance = 0, 0
    lenDoubletDominance, lenSingletDominance = 4000, 4000
    dataDoubletDominance, dataSingletDominance = np.zeros((3, 18, 4000)), np.zeros((3, 18, 4000))
    for angle in [90-10*i for i in range(18)]:
        df_result = calcDoubletFraction(phi, Ca, 1, r, 0, angle, 0)

        # doublet dominance
        if np.mean(df_result[0][0]) > k2:
            lenDoubletDominance = min(df_result[1], lenDoubletDominance) # update union length
            # df
            ax1[0, 0].plot(np.arange(df_result[1])[::slicing], df_result[0][0][:df_result[1]][::slicing], label = "angle = {}".format(angle))
            dataDoubletDominance[0, numberDoubletDominance, :df_result[1]] = df_result[0][0][:df_result[1]]

            # rv
            rv_result = getRelativeViscosity(phi, Ca, 4000, angle, 0)[st:]
            ax1[1, 0].plot(np.arange(len(rv_result))[::slicing], rv_result[::slicing], label = "angle = {}".format(angle))
            dataDoubletDominance[1, numberDoubletDominance, :len(rv_result)] = rv_result

            # iv
            iv_result = getIntrinsicViscosity(phi, Ca, 4000, angle, 0)[st:]
            ax1[2, 0].plot(np.arange(len(iv_result))[::slicing], iv_result[::slicing], label = "angle = {}".format(angle))
            dataDoubletDominance[2, numberDoubletDominance, :len(iv_result)] = iv_result

            numberDoubletDominance += 1

            if set_count == 1:
                angleSelected[set_count].append(angle)

        # singlet dominance
        elif np.mean(df_result[0][0]) < k1:
            lenSingletDominance = min(df_result[1], lenSingletDominance) # update union length
            # df
            ax2[0, 0].plot(np.arange(df_result[1])[::slicing], df_result[0][0][:df_result[1]][::slicing], label = "angle = {}".format(angle))
            dataSingletDominance[0, numberSingletDominance, :df_result[1]] = df_result[0][0][:df_result[1]]

            # rv
            rv_result = getRelativeViscosity(phi, Ca, 4000, angle, 0)[st:]
            ax2[1, 0].plot(np.arange(len(rv_result))[::slicing], rv_result[::slicing], label = "angle = {}".format(angle))
            dataSingletDominance[1, numberSingletDominance, :len(rv_result)] = rv_result

            # iv
            iv_result = getIntrinsicViscosity(phi, Ca, 4000, angle, 0)[st:]
            ax2[2, 0].plot(np.arange(len(iv_result))[::slicing], iv_result[::slicing], label = "angle = {}".format(angle))
            dataSingletDominance[2, numberSingletDominance, :len(iv_result)] = iv_result

            numberSingletDominance += 1

            if set_count != 1:
                angleSelected[set_count].append(angle)

    # error bar plot
    # doublet dominance
    if numberDoubletDominance > 1:
        for i in range(3):
            ax1[i, 1].errorbar(np.arange(lenDoubletDominance)[::slicing],
            np.average(dataDoubletDominance[i, :numberDoubletDominance, :lenDoubletDominance], axis = 0)[::slicing],
            yerr = np.std(dataDoubletDominance[i, :numberDoubletDominance, :lenDoubletDominance], axis = 0)[::slicing])
        
        if set_count == 1:
            timeTrajectories[set_count] = np.average(dataDoubletDominance[:, :numberDoubletDominance, :lenDoubletDominance], axis = 1)
    
    # singlet dominance
    if numberSingletDominance > 1:
        for i in range(3):
            ax2[i, 1].errorbar(np.arange(lenSingletDominance)[::slicing],
            np.average(dataSingletDominance[i, :numberSingletDominance, :lenSingletDominance], axis = 0)[::slicing],
            yerr = np.std(dataSingletDominance[i, :numberSingletDominance, :lenSingletDominance], axis = 0)[::slicing])
        if set_count != 1:
            timeTrajectories[set_count] = np.average(dataSingletDominance[:, :numberSingletDominance, :lenSingletDominance], axis = 1)

    # plot setting
    for i in range(rows):
        for j in range(columns):
            ax1[i, j].set_xlabel("{}".format(r'$\dot \gamma t$'), fontsize = 15)
            ax1[i, j].set_ylabel(symbol_list[i], fontsize = 15)
            ax1[i, j].tick_params(labelsize = 12)
            if numberDoubletDominance > 0 and j == 0:
                ax1[i, j].legend()

            ax2[i, j].set_xlabel("{}".format(r'$\dot \gamma t$'), fontsize = 15)
            ax2[i, j].set_ylabel(symbol_list[i], fontsize = 15)
            ax2[i, j].tick_params(labelsize = 12)
            if numberSingletDominance > 0 and j == 0:
                ax2[i, j].legend()
    fig1.suptitle("phi = {}, Ca = {}: Doublet Dominance\n(k = {}, # of simulations = {})".format(phi, Ca, k2, numberDoubletDominance), fontsize = 20)
    fig1.tight_layout(rect=[0, 0, 1, 0.95])
    fig1.savefig("Pictures/tmp/ThreeTimeSeries_phi_{}_Ca_{}_DoubletDominance_k_{}.png".format(phi, Ca, k2), dpi = 100)
    fig2.suptitle("phi = {}, Ca = {}: Singlet Dominance\n(k = {}, # of simulations = {})".format(phi, Ca, k1, numberSingletDominance), fontsize = 20)
    fig2.tight_layout(rect=[0, 0, 1, 0.95])
    fig2.savefig("Pictures/tmp/ThreeTimeSeries_phi_{}_Ca_{}_SingletDominance_k_{}.png".format(phi, Ca, k1), dpi = 100)

print("Finish finding singlet and doublet dominance")
print("===============================================================================")
print("===============================================================================")
print("Start making three variable plots")

columns, rows = 4, 3

for randomTrial in range(4):
    print("Random Trial {}".format(randomTrial + 1))
    fig, ax = plt.subplots(rows, columns, figsize = (7.5*columns, 6*rows))

    for i in range(columns - 1): # random pick
        for j in range(len(parameter_set)): # parameter set
            [phi, Ca] = parameter_set[j]
            angle = random.choice(angleSelected[j])
            df_result = calcDoubletFraction(phi, Ca, 1, r, 0, angle, 0)
            ax[0, i].plot(np.arange(df_result[1])[::slicing], df_result[0][0][:df_result[1]][::slicing], label = "phi = {}, Ca = {}, angle = {}".format(phi, Ca, angle))
            rv_result = getRelativeViscosity(phi, Ca, 4000, angle, 0)[st:]
            ax[1, i].plot(np.arange(df_result[1])[::slicing], rv_result[:df_result[1]][::slicing], label = "phi = {}, Ca = {}, angle = {}".format(phi, Ca, angle))
            iv_result = getIntrinsicViscosity(phi, Ca, 4000, angle, 0)[st:]
            ax[2, i].plot(np.arange(df_result[1])[::slicing], iv_result[:df_result[1]][::slicing], label = "phi = {}, Ca = {}, angle = {}".format(phi, Ca, angle))

    # ensemble averaged
    for i in range(3): # variable
        for j in range(len(parameter_set)): # parameter set
            [phi, Ca] = parameter_set[j]
            ts = timeTrajectories[j][i, :]
            ax[i, columns-1].plot(np.arange(len(ts))[::slicing], ts[::slicing], label = "phi = {}, Ca = {}".format(phi, Ca))

    # plot setting
    for j in range(columns):
        ax[0, j].set_title("random {}".format(j+1) if j != (columns - 1) else "Ensemble Averaged", fontsize = 20)
        for i in range(rows):
            ax[i, j].set_xlabel("{}".format(r'$\dot \gamma t$'), fontsize = 15)
            ax[i, j].set_ylabel(symbol_list[i], fontsize = 15)
            ax[i, j].tick_params(labelsize = 12)
            ax[i, j].legend()

    fig.savefig("Pictures/tmp/ThreeTimeSeries_Comparison_k1_{}_k2_{}_random_{}.png".format(k1, k2, randomTrial+1), dpi = 100)
print("Finish making three variable plots")
print("===============================================================================")
print('Total time elapsed = {}'.format(str(datetime.timedelta(seconds=time.time()-start_time))))
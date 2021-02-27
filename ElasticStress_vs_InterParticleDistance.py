# ===============================================================================
# Copyright 2021 An-Jun Liu
# Last Modified Date: 01/30/2021
# ===============================================================================
import numpy as np 
import matplotlib.pyplot as plt
from RBC_Utilities import calcDoubletFraction, getIntrinsicViscosity, getStress
import scipy as scipy

"""
This code is to plot elastic stress tensor vs interparticle distance, avg or max, in doublet state.
Also provided linear regression line.
"""

def plot():
    result = calcDoubletFraction(phi, Ca, 1.0, r, 0, angle, 0)
    dfg = result[0][0]
    diffpos = result[3][0,0:len(dfg)]
    stress = getStress(phi, Ca, 0, 0, ncycle, angle, 0)[:,3][0:len(dfg)]
    time = np.array(list(range(len(dfg))))

    w = 3
    st, et = 0, 1900
    local_min_time = []

    avg_distance = []
    max_distance = []
    avg_stress = []
    max_stress = []

    num_min = 0
    threshold_distance = 0.5
    threshold_stress = 0.0001
    distance_sum, distance_period, distance_max = 0, 0, -1
    stress_sum, stress_max = 0, -1

    min_threshold = 1
    num_min_stress, stress_min_num = 0, 0

    for t in time[st+w:et-w]:
        if dfg[t] == 1: # in doubelt state
            # interval
            if num_min >= min_threshold: # start storing values
                distance_period += 1
                # avg
                distance_sum += diffpos[t]
                stress_sum += stress[t]
                # max
                if diffpos[t] > distance_max:
                    distance_max = diffpos[t]
                if stress[t] > stress_max:
                    stress_max = stress[t]


            # local min
            if (((diffpos[t-1]-diffpos[t]) > 0) and ((diffpos[t+1]-diffpos[t]) > 0)):
                if (((diffpos[t-w]-diffpos[t]) > threshold_distance) and ((diffpos[t+w]-diffpos[t]) > threshold_distance)):
                    if num_min < min_threshold:
                        num_min += 1
                    else:
                        avg_distance.append(distance_sum/distance_period)
                        max_distance.append(distance_max)
                        avg_stress.append(stress_sum/distance_period)
                        max_stress.append(stress_max)
        
                        distance_sum, distance_period, distance_max = 0, 0, -1
                        stress_sum, stress_max = 0, -1

            if (((stress[t-1]-stress[t])>0) and ((stress[t+1]-stress[t])>0)):
                if (((stress[t-w]-stress[t])>threshold_stress) and ((stress[t+w]-stress[t])>threshold_stress)):
                    if num_min_stress < min_threshold:
                        num_min_stress += 1
                    else:
                        stress_min_num += 1

        else:
            num_min = 0
            num_min_stress = 0

    slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(avg_distance, avg_stress)
    fit_dis = np.array(avg_distance)
    fit_stress = fit_dis*slope + intercept

    plt.figure(figsize=(6, 6))
    plt.title('phi = {}, Ca = {}, angle = {}, criteria r = {}Dm\n'.format(phi, Ca, angle, r), fontsize = 15)
    plt.plot(avg_distance, avg_stress, 'o', label = "Data")
    plt.plot(fit_dis, fit_stress, label = "Fitted line")
    plt.legend()
    #plt.text(8.5, 0.00047, "y={}x+{}\n".format(np.round(slope,6), np.round(intercept,6)) + r"$R^{2}$" + "={}".format(np.round(r_value**2,6)))
    plt.annotate("y={}x+{}\n".format(np.round(slope,7), np.round(intercept,6)) + r"$R^{2}$" + "={}".format(np.round(r_value**2,6)), xy=(0.55, 0.75), xycoords='axes fraction')
    plt.xlabel("interparticle distance", fontsize = 15)
    plt.ylabel(r'$\sigma _{elas, yx}$', fontsize = 15)
    plt.tight_layout()
    #plt.savefig("./Pictures/demo5.png")
    plt.show()

ncycle = 2000
r = 1
phi, Ca, angle = 3.8, 0.02, 40
plot()
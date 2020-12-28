import numpy as np 
import matplotlib.pyplot as plt
from Doublet_Functions import calcDoubletFraction, getInstrinsicViscosity, getStress
import scipy as scipy

"""
This code is to plot elastic stress tensor vs interparticle distance, avg or max.
Also provided linear regression line.
"""

def plot(plot_or_not):
    w = 3
    st, et = 0, 1900

    avg_distance = []
    avg_stress = []

    for i in range(18):
        angle = 90 - 10*i

        result = calcDoubletFraction(phi, Ca, 1.0, r, 0, angle)
        dfg = result[0][0]
        diffpos = result[3][0,0:len(dfg)]
        stress = getStress(phi, Ca, 0, 0, ncycle, angle)[:,3][0:len(dfg)]
        time = np.array(list(range(len(dfg))))

        num_min = 0
        distance_sum, distance_period = 0, 0
        stress_sum = 0

        min_threshold = 3

        for t in time[st+w:et-w]:
            if dfg[t] == 1: # in doubelt state
                # interval
                if num_min >= min_threshold: # start storing values
                    distance_period += 1
                    # avg
                    distance_sum += diffpos[t]
                    stress_sum += stress[t]


                # local min
                if (((diffpos[t-1]-diffpos[t]) > 0) and ((diffpos[t+1]-diffpos[t]) > 0)):
                    if (((diffpos[t-w]-diffpos[t]) > 0) and ((diffpos[t+w]-diffpos[t]) > 0)):
                        if num_min < min_threshold:
                            num_min += 1
                        else:
                            avg_distance.append(distance_sum/distance_period)
                            avg_stress.append(stress_sum/distance_period)
            
                            distance_sum, distance_period = 0, 0
                            stress_sum = 0

            else:
                num_min = 0

    slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(avg_distance, avg_stress)
    fit_dis = np.array(avg_distance)
    fit_stress = fit_dis*slope + intercept

    plt.figure(figsize=(6, 6))
    plt.title('phi = {}, Ca = {}, criteria r = {}Dm\n'.format(phi, Ca, r), fontsize = 15)
    plt.plot(avg_distance, avg_stress, 'o', label = "Data")
    plt.plot(fit_dis, fit_stress, label = "Fitted line")
    plt.legend()
    plt.annotate("y={}x+{}\n".format(np.round(slope,7), np.round(intercept,6)) + r"$R^{2}$" + "={}".format(np.round(r_value**2,6)), xy=(0.55, 0.75), xycoords='axes fraction')
    plt.xlabel("interparticle distance", fontsize = 15)
    plt.ylabel(r'$\sigma _{elas, yx}$', fontsize = 15)
    plt.tight_layout()
    if plot_or_not == 1:
        plt.savefig("/Users/andrewliu/Code/YL_project/Pictures/ElasticStress_yx_vs_interparticledistance/phi_{}_Ca_{}_r_{}.png".format(phi, Ca, r))
    else:
        plt.show()

ncycle = 2000
r = 1
#Phi, Ca, angle = 6.0, 0.08, -70
phi, Ca = 6.0, 0.06
plot(1)
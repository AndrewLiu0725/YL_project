# ===============================================================================
# Copyright 2021 An-Jun Liu
# Last Modified Date: 03/03/2021
# ===============================================================================
import numpy as np
import matplotlib.pyplot as plt
import time
import datetime
from RBC_Utilities import calcDoubletFraction, getIntrinsicViscosity
import os
import sys

"""
This code is to compare the instability of the two-cell system before and after increasing the simulation time to 4000 strains.
A plot containing a histogram and a cmf is provided.
"""

start_time = time.time()

root_folder = "/userdata4/ajliu/RBC_doublet/"

ncycle = 4000
vol = 746.3163
phi_range =[]
for x in range(30, 41):
    phi = 2*vol/(24*x**2)
    phi_range.append(float(str(phi*100)[:3]))
Ca_range = [i*0.01 for i in range(1, 21)]
angle_range = [90 - 10*j for j in range(18)]

deviation_origin = []
deviation_extend = []

for phi in phi_range:
    for Ca in Ca_range:
        for angle in angle_range:
            job = "phi_{}_Re_0.1_Ca_{}_aggregation_1KT_ncycle_{}_np_2_angle_{}".format(phi, Ca, ncycle, angle)
            if os.path.isfile(root_folder+"Data/"+job+"/data/bond0_t{}.vtk".format(int(3669*2*(int(ncycle/2)-1)))):
                try:
                    iv = getIntrinsicViscosity(phi, Ca, ncycle, angle, 0)
                    first_extend = np.mean(iv[int(0.5*ncycle):int(0.75*ncycle)])
                    second_extend = np.mean(iv[int(0.75*ncycle):ncycle]) 

                    first_origin = np.mean(iv[int(0.5*int(ncycle/2)):int(0.75*int(ncycle/2))])
                    second_origin = np.mean(iv[int(0.75*int(ncycle/2)):int(ncycle/2)]) 

                    deviation_extend.append(abs((second_extend-first_extend)/first_extend))
                    deviation_origin.append(abs((second_origin-first_origin)/first_origin))

                except KeyboardInterrupt:
                    print("Pressed ctrl+C\n")
                    sys.exit()

                except Exception as e:
                    print("Error: Case (phi = {}, Ca = {}, angle = {}): \n{}".format(phi, Ca, angle, e))


data_num = len(deviation_origin)
bin_num = 10
upper = max(max(deviation_extend), max(deviation_origin))
bin_width = (upper-0)/bin_num
edge = np.linspace(0, upper, num = (bin_num+1))
center = np.array([(edge[i]+edge[i+1])/2 for i in range(bin_num)])
density_origin = np.zeros(bin_num)
density_extend = np.zeros(bin_num)

deviation_extend.sort()
deviation_origin.sort()

ptr_origin = 0
for dev in deviation_origin:
    # find ptr such that edge[ptr_orgin] < dev <= edge[ptr_orgin+1]
    while dev > edge[ptr_origin+1]:
        ptr_origin += 1
    density_origin[ptr_origin] += 1 

ptr_extend = 0
for dev in deviation_extend:
    # find ptr such that edge[ptr_orgin] < dev <= edge[ptr_orgin+1]
    while dev > edge[ptr_extend+1]:
        ptr_extend += 1
    density_extend[ptr_extend] += 1 

density_extend = density_extend/data_num
density_origin = density_origin/data_num

cmf_origin = np.zeros(bin_num)
cmf_origin[0] = density_origin[0]
for i in range(1, bin_num):
    cmf_origin[i] = cmf_origin[i-1] + density_origin[i]

cmf_extend = np.zeros(bin_num)
cmf_extend[0] = density_extend[0]
for i in range(1, bin_num):
    cmf_extend[i] = cmf_extend[i-1] + density_extend[i]

fig, (ax1, ax2) = plt.subplots(1, 2)
bar_width = 3/8
fig.suptitle("Comparison of instability for two-cell system \nafter increasing simulation time\n")
ax1.bar(center-bin_width*bar_width/2, density_origin, width = bin_width*bar_width, label = "2000 strains", color = 'r')
ax1.bar(center+bin_width*bar_width/2, density_extend, width = bin_width*bar_width, label = "4000 strains", color = 'b')
ax1.set_title("Histogram")
ax1.set(xlabel = "Instability", ylabel = "Density")
ax1.legend()

ax2.plot(center, cmf_origin, label = "2000 strains", color = 'r')
ax2.plot(center, cmf_extend, label = "4000 strains", color = 'b')
ax2.set_title("Cumulative mass function")
ax2.set(xlabel = "Instability", ylabel = "Density")
ax2.legend()

fig.tight_layout()
plt.subplots_adjust(top = 0.85)
plt.savefig("./Pictures/TwoCellSystem_Instability_Histogram_Comparison.png", dpi = 300)

print('Total time elapsed = {}'.format(str(datetime.timedelta(seconds=time.time()-start_time))))
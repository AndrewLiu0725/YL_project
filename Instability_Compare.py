# ===============================================================================
# Copyright 2021 An-Jun Liu
# Last Modified Date: 03/13/2021
# ===============================================================================
import numpy as np
import matplotlib.pyplot as plt
import time
import datetime
from RBC_Utilities import calcDoubletFraction, getIntrinsicViscosity, getSuspensionParameterSets
import os
import sys

"""
This code is to compare the instability of the two-cell system before and after increasing the simulation time to 4000 strains.
Added the data of suspension system to see the difference.
A plot containing a histogram and a cmf is provided.
"""

start_time = time.time()


## two-cell system
# ===============================================================================
print("Start colleting the data for two-cell system.")
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
deviation_extend_fix_500 = []
deviation_extend_fix_ratio = []

for phi in phi_range:
    for Ca in Ca_range:
        for angle in angle_range:
            job = "phi_{}_Re_0.1_Ca_{}_aggregation_1KT_ncycle_{}_np_2_angle_{}".format(phi, Ca, ncycle, angle)
            if os.path.isfile(root_folder+"Data/"+job+"/data/bond0_t{}.vtk".format(int(3669*2*(int(ncycle/2)-1)))):
                try:
                    iv = getIntrinsicViscosity(phi, Ca, ncycle, angle, 0)
                    first_extend_fix_ratio = np.mean(iv[int(0.5*ncycle):int(0.75*ncycle)])
                    second_extend_fix_ratio = np.mean(iv[int(0.75*ncycle):ncycle]) 

                    first_extend_fix_500 = np.mean(iv[-1000:-500])
                    second_extend_fix_500 = np.mean(iv[-500:])

                    first_origin = np.mean(iv[int(0.5*int(ncycle/2)):int(0.75*int(ncycle/2))])
                    second_origin = np.mean(iv[int(0.75*int(ncycle/2)):int(ncycle/2)]) 

                    deviation_extend_fix_ratio.append(abs((second_extend_fix_ratio-first_extend_fix_ratio)/first_extend_fix_ratio))
                    deviation_extend_fix_500.append(abs((second_extend_fix_500-first_extend_fix_500)/first_extend_fix_500))
                    deviation_origin.append(abs((second_origin-first_origin)/first_origin))

                except KeyboardInterrupt:
                    print("Pressed ctrl+C\n")
                    sys.exit()

                except Exception as e:
                    print("Error: Case (phi = {}, Ca = {}, angle = {}): \n{}".format(phi, Ca, angle, e))



## suspension system
# ===============================================================================
print("Start colleting the data for suspension system.")

[phis, parameter_set] = getSuspensionParameterSets()

# calculate deviation
deviation_suspension = []
max_timesteps = 10000

for phi in phis:
    for Ca in parameter_set[phi].keys():
        for ensemble_id in parameter_set[phi][Ca]:
            try:
                iv = getIntrinsicViscosity(phi, Ca, max_timesteps, ensemble_id, 1)
                first = np.mean(iv[int(0.5*len(iv)):int(0.75*len(iv))])
                second = np.mean(iv[int(0.75*len(iv)):len(iv)]) 
                deviation_suspension.append(abs((second-first)/first))

            except KeyboardInterrupt:
                print("Pressed ctrl+C\n")
                sys.exit()

            except Exception as e:
                print("Error: Case (phi = {}, Ca = {}, ensemble_id = {}): \n{}".format(phi, Ca, ensemble_id, e))



## make plot
# ===============================================================================
bin_num = 10
upper = max(max(deviation_extend_fix_500), max(deviation_extend_fix_ratio), max(deviation_origin), max(deviation_suspension))
bin_width = (upper-0)/bin_num
edge = np.linspace(0, upper, num = (bin_num+1))
center = np.array([(edge[i]+edge[i+1])/2 for i in range(bin_num)])
density_origin = np.zeros(bin_num)
density_extend_fix_500 = np.zeros(bin_num)
density_extend_fix_ratio = np.zeros(bin_num)
density_suspension = np.zeros(bin_num)

deviation_extend_fix_ratio.sort()
deviation_extend_fix_500.sort()
deviation_origin.sort()
deviation_suspension.sort()

ptr_origin = 0
for dev in deviation_origin:
    # find ptr such that edge[ptr_orgin] < dev <= edge[ptr_orgin+1]
    while dev > edge[ptr_origin+1]:
        ptr_origin += 1
    density_origin[ptr_origin] += 1 

ptr_extend_fix_500 = 0
for dev in deviation_extend_fix_500:
    # find ptr such that edge[ptr_orgin] < dev <= edge[ptr_orgin+1]
    while dev > edge[ptr_extend_fix_500+1]:
        ptr_extend_fix_500 += 1
    density_extend_fix_500[ptr_extend_fix_500] += 1 

ptr_extend_fix_ratio = 0
for dev in deviation_extend_fix_ratio:
    # find ptr such that edge[ptr_orgin] < dev <= edge[ptr_orgin+1]
    while dev > edge[ptr_extend_fix_ratio+1]:
        ptr_extend_fix_ratio += 1
    density_extend_fix_ratio[ptr_extend_fix_ratio] += 1 

ptr_suspension = 0
for dev in deviation_suspension:
    # find ptr such that edge[ptr_orgin] < dev <= edge[ptr_orgin+1]
    while dev > edge[ptr_suspension+1]:
        ptr_suspension += 1
    density_suspension[ptr_suspension] += 1 

density_extend_fix_500 = density_extend_fix_500/len(deviation_origin)
density_extend_fix_ratio = density_extend_fix_ratio/len(deviation_origin)
density_origin = density_origin/len(deviation_origin)
density_suspension = density_suspension/len(deviation_suspension)

cmf_origin = np.zeros(bin_num)
cmf_origin[0] = density_origin[0]
for i in range(1, bin_num):
    cmf_origin[i] = cmf_origin[i-1] + density_origin[i]

cmf_extend_fix_500 = np.zeros(bin_num)
cmf_extend_fix_500[0] = density_extend_fix_500[0]
for i in range(1, bin_num):
    cmf_extend_fix_500[i] = cmf_extend_fix_500[i-1] + density_extend_fix_500[i]

cmf_extend_fix_ratio = np.zeros(bin_num)
cmf_extend_fix_ratio[0] = density_extend_fix_ratio[0]
for i in range(1, bin_num):
    cmf_extend_fix_ratio[i] = cmf_extend_fix_ratio[i-1] + density_extend_fix_ratio[i]

cmf_suspension = np.zeros(bin_num)
cmf_suspension[0] = density_suspension[0]
for i in range(1, bin_num):
    cmf_suspension[i] = cmf_suspension[i-1] + density_suspension[i]

fig, (ax1, ax2) = plt.subplots(1, 2)
bar_width = 1/4
fig.suptitle("Comparison of instability for two-cell system \nafter increasing simulation time\n")
ax1.bar(center - (3/2)*bin_width*bar_width, density_origin, width = bin_width*bar_width, label = "2000 strains", color = 'r')
ax1.bar(center - (1/2)*bin_width*bar_width, density_extend_fix_ratio, width = bin_width*bar_width, label = "4000 strains (quarter)", color = 'b')
ax1.bar(center + (1/2)*bin_width*bar_width, density_extend_fix_500, width = bin_width*bar_width, label = "4000 strains (500 strains)", color = 'k')
ax1.bar(center + (3/2)*bin_width*bar_width, density_suspension, width = bin_width*bar_width, label = "suspension", color = 'g')
ax1.set_title("Histogram")
ax1.set(xlabel = "Instability", ylabel = "Density")
ax1.legend()

ax2.plot(center, cmf_origin, label = "2000 strains", color = 'r')
ax2.plot(center, cmf_extend_fix_ratio, label = "4000 strains (quarter)", color = 'b')
ax2.plot(center, cmf_extend_fix_500, label = "4000 strains (500 strains)", color = 'k')
ax2.plot(center, cmf_suspension, label = "suspension", color = 'g')
ax2.set_title("Cumulative mass function")
ax2.set(xlabel = "Instability", ylabel = "Density")
ax2.legend()

fig.tight_layout()
plt.subplots_adjust(top = 0.85)
plt.savefig("./Pictures/TwoCellSystem_Instability_Histogram_Comparison.png", dpi = 300)

print('Total time elapsed = {}'.format(str(datetime.timedelta(seconds=time.time()-start_time))))
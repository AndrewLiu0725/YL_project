import numpy as np 
import matplotlib.pyplot as plt
from Doublet_Functions import calcDoubletFraction, getInstrinsicViscosity, getRelativeViscosity
import time as time
import datetime as datetime

"""
This code is to plot intrinsic and relative viscosity vs Ca 
for both doublet and two-singlet state with all volume fraction.
"""

### Calculation here

def clac_stat(time_window, ratio_k):
    global stat
    start_timestep = 10
    result = calcDoubletFraction(phi, Ca, 1.0, r, 0, angle)
    state = result[2][0][0, :2000] # end at 2000 if the total number of timesteps exceeds 2000
    end_timestep = len(state)
    iv = getInstrinsicViscosity(phi, Ca, 0, ncycle, angle)[:end_timestep]
    rv = getRelativeViscosity(phi, Ca, 0, ncycle, angle)[:end_timestep]
    target_count = time_window*ratio_k
    
    # calculate the stat
    # initialization the count
    count = np.zeros(4)
    half_width = int(time_window/2)
    t = start_timestep + half_width
    for i in range(start_timestep, start_timestep + time_window):
        count[state[i]-1] += 1
    # run over both states
    for j in range(2):
        if (count[j] > target_count):
            stat[0, j, 0] += iv[t]
            stat[1, j, 0] += (iv[t])**2
            stat[2, j, 0] += 1
            stat[0, j, 1] += rv[t]
            stat[1, j, 1] += (rv[t])**2
            stat[2, j, 1] += 1
    
    for t in range(start_timestep + half_width + 1, end_timestep - half_width):
        count[state[t - half_width - 1] - 1] -= 1
        count[state[t + half_width] - 1] += 1
        # run over both states
        for j in range(2):
            if (count[j] > target_count):
                stat[0, j, 0] += iv[t]
                stat[1, j, 0] += (iv[t])**2
                stat[2, j, 0] += 1
                stat[0, j, 1] += rv[t]
                stat[1, j, 1] += (rv[t])**2
                stat[2, j, 1] += 1



start_time = time.time()

# Parameters
r = 1.0
ncycle = 2000
time_window = 100
ratio_k = 0.9
min_count = 100

vol = 746.3163
phis  =[]
for x in range(30, 41):
    phi = 2*vol/(24*x**2)
    phis.append(float(str(phi*100)[:3]))
phi_range = np.array(phis)
Ca_range = np.array([i*0.01 for i in range(1, 21)])

# (iv/rv, d/s, mean/std, phi, Ca)
data = np.zeros((2, 2, 2, len(phi_range), len(Ca_range)))

for phi_index, phi in enumerate(phi_range):
    print(phi)
    for Ca_index, Ca in enumerate(Ca_range):
        ## Store stat for each set of phi and Ca: x/x^2/count, d/k, iv/rv
        stat = np.zeros((3, 2, 3))
        for j in range(18):
            angle = 90-10*j
            clac_stat(time_window, ratio_k)

        ## store the value       
        for i in range(2): # run over state category
            if  stat[2, i, 0] > min_count: # this state exists
                data[0, i, 0, phi_index, Ca_index] = stat[0, i, 0]/stat[2, i, 0]
                data[0, i, 1, phi_index, Ca_index] = np.sqrt(stat[1, i, 0]/stat[2, i, 0] - (stat[0, i, 0]/stat[2, i, 0])**2)
                data[1, i, 0, phi_index, Ca_index] = stat[0, i, 1]/stat[2, i, 1]
                data[1, i, 1, phi_index, Ca_index] = np.sqrt(stat[1, i, 1]/stat[2, i, 1] - (stat[0, i, 1]/stat[2, i, 1])**2)


### Plot here
np.save("IRstress_data.npy", data)

### vs Ca
import matplotlib.colors as mcolors
import random
cmap = random.sample(list(mcolors.CSS4_COLORS.keys()), len(phi_range))

for i in range(2): # run over iv and rv
    plt.figure(figsize=(12,9))
    #cmap = ['b','g','r','c','m','y','k']
    for j in range(2): # run over doublet and kayaking state
        for phi_index, phi in enumerate(phi_range):
            nz_range = np.nonzero(data[i, j, 0, phi_index, :])
            if j == 0:
                plt.plot(Ca_range[nz_range], data[i, j, 0, phi_index, :][nz_range], label = 'phi = {}, doublet'.format(phi), color = cmap[phi_index])
            elif j == 1:       
                plt.plot(Ca_range[nz_range], data[i, j, 0, phi_index, :][nz_range], label = 'phi = {}, singlet'.format(phi), ls = '--', color = cmap[phi_index])    
            
    plt.xlabel('Ca', fontsize = 20)
    if i == 0:
        text = "Intrinsic Viscosity"
    elif i == 1:
        text = "Relative viscosity"
    plt.ylabel(text, fontsize = 20)
    plt.title(text+' vs Ca', fontsize=30)
    plt.legend()
    plt.savefig('./Pictures/TwoCellSystem_{}_Ca_summary_r_{}_timewindow_{}_k_{}.png'.format("".join(text.split()), r, time_window, ratio_k), dpi = 300)


### vs phi
for Ca_index, Ca in enumerate(Ca_range):
    fig, ax1 = plt.subplots(figsize = (9, 9))
    color = 'tab:red'
    ax1.set_xlabel(r'$\phi$', fontsize = 20)
    ax1.set_ylabel("Intrinsic Viscosity", color = color, fontsize = 20)
    ax1.set_ylim([min(np.min(data[0, 0, 0, :, :][np.nonzero(data[0, 0, 0, :, :])]), np.min(data[0, 1, 0, :, :][np.nonzero(data[0, 1, 0, :, :])]))*0.9,
     max(np.max(data[0, 0, 0, :, :][np.nonzero(data[0, 0, 0, :, :])]), np.max(data[0, 1, 0, :, :][np.nonzero(data[0, 1, 0, :, :])]))*1.1])
    ax1.set_title("Ca = {}".format(Ca), fontsize = 30)
    ax1.tick_params(axis='y', labelcolor = color, labelsize = 20)
    ax1.tick_params(axis='x', labelsize = 20)
    
    for i in range(2): # run over state
        nz_range = np.nonzero(data[0, i, 0, :, Ca_index])
        if len(nz_range) > 0:
            if i == 0:
                p1, = ax1.plot(0.01*phi_range[nz_range], data[0, i, 0, :, Ca_index][nz_range], color = color, label = "Intrinsic Viscosity, doublet")
            else:
                p2, = ax1.plot(0.01*phi_range[nz_range], data[0, i, 0, :, Ca_index][nz_range], color = color, ls = '--', label = "Intrinsic Viscosity, singlet")  

    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

    color = 'tab:blue'
    ax2.set_ylabel("Relative Viscosity", color = color, fontsize = 20)
    ax2.set_ylim([min(np.min(data[1, 0, 0, :, :][np.nonzero(data[1, 0, 0, :, :])]), np.min(data[1, 1, 0, :, :][np.nonzero(data[1, 1, 0, :, :])]))*0.9,
     max(np.max(data[1, 0, 0, :, :][np.nonzero(data[1, 0, 0, :, :])]), np.max(data[1, 1, 0, :, :][np.nonzero(data[1, 1, 0, :, :])]))*1.1])
    for i in range(2): # run over state
        nz_range = np.nonzero(data[1, i, 0, :, Ca_index])
        if len(nz_range) > 0:
            if i == 0:
                p3, = ax2.plot(0.01*phi_range[nz_range], data[1, i, 0, :, Ca_index][nz_range], color = color, label = "Intrinsic Viscosity, doublet")
            else:
                p4, = ax2.plot(0.01*phi_range[nz_range], data[1, i, 0, :, Ca_index][nz_range], color = color, ls = '--', label = "Intrinsic Viscosity, singlet")  

    ax2.tick_params(axis = 'y', labelcolor = color, labelsize = 20)

    lines = [p1, p2, p3, p4]
    ax1.legend(lines, [l.get_label() for l in lines], bbox_to_anchor=(1.0,1.2))
    fig.tight_layout()  # otherwise the right y-label is slightly clipped
    fig.savefig('./Pictures/TwoCellSystem_Viscosity_vs_phi_Ca_{}_fixed_ylim_r_{}_timewindow_{}_k_{}.png'.format(Ca, r, time_window, ratio_k), dpi = 100)
    ax2.set_ylim([min(np.min(data[1, 0, 0, :, Ca_index][np.nonzero(data[1, 0, 0, :, Ca_index])]), np.min(data[1, 1, 0, :, Ca_index][np.nonzero(data[1, 1, 0, :, Ca_index])])),
     max(np.max(data[1, 0, 0, :, Ca_index][np.nonzero(data[1, 0, 0, :, Ca_index])]), np.max(data[1, 1, 0, :, Ca_index][np.nonzero(data[1, 1, 0, :, Ca_index])]))])
    ax1.set_ylim([min(np.min(data[0, 0, 0, :, Ca_index][np.nonzero(data[0, 0, 0, :, Ca_index])]), np.min(data[0, 1, 0, :, Ca_index][np.nonzero(data[0, 1, 0, :, Ca_index])])),
     max(np.max(data[0, 0, 0, :, Ca_index][np.nonzero(data[0, 0, 0, :, Ca_index])]), np.max(data[0, 1, 0, :, Ca_index][np.nonzero(data[0, 1, 0, :, Ca_index])]))])
    fig.savefig('./Pictures/TwoCellSystem_Viscosity_vs_phi_Ca_{}_r_{}_timewindow_{}_k_{}.png'.format(Ca, r, time_window, ratio_k), dpi = 300)

print('Total time elapsed = {}'.format(str(datetime.timedelta(seconds=time.time()-start_time))))
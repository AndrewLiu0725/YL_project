# ===============================================================================
# Copyright 2021 An-Jun Liu
# Last Modified Date: 08/06/2022
# ===============================================================================
import numpy as np
from scipy import stats
import math
import matplotlib.pyplot as plt 
from matplotlib.ticker import MultipleLocator
from RBC_Utilities import getSuspensionParameterSets, calcDoubletFraction
import time
import datetime
import pickle

"""
This code is to make the system's roation time (seperate doublet and siglets state) vs Ca 
with varied volume fractions for suspension and two-cell system.
"""

# setting
# ===============================================================================
path_preprocess = "/userdata4/ajliu/Data_Transfer/"
path_CT = "/raid6/ctliao/Data/HI_ordering/"
path_AJ = "/userdata4/ajliu/RBC_doublet/Data/"

SAVE = 0
READ = 1
SHOW = 0
print("Setting:\nSAVE:{}, READ:{}, SHOW:{}".format(SAVE, READ, SHOW))
window_width = 100
threshold = 0.9
scale_factor = [1, 4000/3669]


# main function
# ===============================================================================
def getRotationTime(input_phi, input_Ca, input_criteria_T, input_criteria_Dms, depend, system):
    [doublet_or_not, indice_pairs, number_of_pairs, particle_numbers, timesteps, Ypos_t, cutoff_frequency, points_per_particle, interval] = calcDoubletFraction(input_phi, input_Ca, input_criteria_T, input_criteria_Dms, 0, depend, system, outputDataType=2)
          
    state = np.zeros((particle_numbers, timesteps))
    target1, target2 = int(window_width*threshold), int(window_width*(1-threshold))
    scale = timesteps/interval

    # mark the doublet state for each particle time series
    for i in range(number_of_pairs):
        j, k = indice_pairs[i, 0], indice_pairs[i, 1]

        # go through the i-th pair
        for t in range(timesteps):
            if doublet_or_not[i, t] > 0:
                state[j, t] = 1
                state[k, t] = 1

    # compute the rotation times of doublet and singlet state dominant time intervals for each particle
    t_rot_stat = [[], []] # [[singlet], [doublet]]

    for i in range(particle_numbers):
        # split the time series according to the state
        # initialize
        split = []
        record = 0 # 1 means in dominance interval, so 0->1->0 specifies a dominace interval
        current_state = 0 # 1 for doublet dominat interval, 0 for singlet one

        dc = np.sum(state[i, :window_width]) # doublet data point count
        if dc > target1:
            record = 1
            current_state = 1
            st = 0
        elif dc < target2:
            record = 1
            current_state = 0
            st = 0

        # go through the time series:
        for t in range(window_width, timesteps-window_width):
            dc += (state[i, t] - state[i, t-window_width]) # update dc
            if record:
                # 1->0
                if (current_state == 1 and dc < target1) or (current_state == 0 and dc > target2):
                    if (t - st) >= window_width:
                        split.append([st, t, current_state])
                    record = 0
            else:
                # 0->1
                if dc >= target1:
                    st = t
                    current_state = 1
                    record = 1
                elif dc <= target2:
                    st = t
                    current_state = 0
                    record = 1

        # edge case (last interval)
        if record and (t - st) >= window_width:
            split.append([st, t, current_state])


        # calculate rotation time
        for slice in split:
            for j in range(points_per_particle):
                node_id = i*points_per_particle + j
                cutoff = int(cutoff_frequency*(slice[1]-slice[0]) + 1)
                P = np.fft.rfft((Ypos_t[int(slice[0]/scale):int(slice[1]/scale), node_id] - np.mean(Ypos_t[int(slice[0]/scale):int(slice[1]/scale), node_id]))) # remove the DC term
                t_rot_stat[slice[2]].append(scale_factor[system]*(slice[1]-slice[0])/(np.argmax(np.abs(P[cutoff:]))+cutoff))
    
    exclude_ratio = 0.2
    # [singlet, doublet]
    return [stats.trim_mean(t_rot_stat[0], exclude_ratio), stats.trim_mean(t_rot_stat[1], exclude_ratio)] 



### main code
# ===============================================================================
start_time = time.time()
if SAVE:
    data = {"TwoCell":{}, "Suspension":{}}

if READ:
    with open("Data/T_rot_TCandS_w_{}_k_{}.pickle".format(window_width, threshold), 'rb') as handle:
        data = pickle.load(handle) # data[system][phi] = [Cas, t_rot, t_rot_std]

criteria_1 = 3
color_list = ['orange', 'green']

plt.rcParams['font.family'] = 'DeJavu Serif'
plt.rcParams['font.serif'] = ['Times New Roman']
# initialize the plots

fig, ax = plt.subplots(figsize = (9, 6))

# run over two-cell and suspension system
# ===============================================================================
for system in ['TwoCell', 'Suspension']:
    # set up the parameter set
    # two-cell
    if system == 'TwoCell':
        vol = 746.3163
        phi_range =[]
        true_phi_range = []
        for x in range(30, 41):
            phi = 2*vol/(24*x**2)
            phi_range.append(float(str(phi*100)[:3]))
            true_phi_range.append(round(phi*100, 1))
        Ca_range = [i*0.01 for i in range(1, 21)]
        angle_range = [90 - 10*j for j in range(18)]
        phi_range.sort()
        true_phi_range.sort()
        chosen_phi = [4.0, 6.0]
    # suspension
    else:
        if SAVE:
            [phi_range, parameter_set] = getSuspensionParameterSets() # need to know the ensemble id

        phi_range = [2.3994, 2.9993, 3.4492, 3.8991, 4.9488, 5.9986]
        Ca_range = [0.01, 0.02, 0.03, 0.06, 0.07, 0.08, 0.09, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2]
        chosen_phi = [3.8991, 5.9986]
        true_phi_range = phi_range

    # run over phi
    for phi_index, phi in enumerate(phi_range):
        if not phi in chosen_phi: continue

        # collect data for this phi
        if READ:
            t_rot = data[system][phi][1]
            t_rot_std = data[system][phi][2]
            Cas = data[system][phi][0]

        else:
            t_rot = [[], []]
            t_rot_std = [[], []]
            Cas = [[], []]

            for Ca_index, Ca in enumerate(Ca_range):
                tmp_t_rot = [[], []]

                # run over ensembles
                # Two-cell
                if system == 'TwoCell':
                    for angle_index, angle in enumerate(angle_range):
                        result = getRotationTime(phi, Ca, 1, 1, angle, 0)
                        for i in range(2):
                            if not math.isnan(result[i]):
                                tmp_t_rot[i].append(result[i])
                # Suspension
                else:
                    for ensemble_id in parameter_set[phi][Ca]:
                        try:
                            result = getRotationTime(phi, Ca, 1, 1, ensemble_id, 1)
                            for i in range(2):
                                if not math.isnan(result[i]):
                                    tmp_t_rot[i].append(result[i])
                        except:
                            pass


                for i in range(2):
                    if len(tmp_t_rot[i]) > criteria_1:
                        t_rot[i].append(np.mean(tmp_t_rot[i]))
                        t_rot_std[i].append(np.std(tmp_t_rot[i]))
                        Cas[i].append(Ca)

        if SAVE:
            data[system][phi] = [Cas, t_rot, t_rot_std] # data[system][phi][Cas/t_rot/t_rot_std][singlet/doublet]

        # plot for this phi
        for state in range(2):
            ax.scatter(Cas[state], t_rot[state], marker='o'if system=='Suspension' else 's',
            facecolors=color_list[chosen_phi.index(phi)]if state else 'none', edgecolors=color_list[chosen_phi.index(phi)],
            label='{}={:.2}%-{}'.format(r'$\phi$', true_phi_range[phi_index], "doublet"if state else 'singlet'))
            #label='{} System\n{}={:.2}%-{}'.format('Two-Cell'if system=='TwoCell' else 'Suspension', r'$\phi$', true_phi_range[phi_index], "doublet"if state else 'singlet')


for i in range(2): # run over singlet and doublet
    y, x = [], []
    for system in ['TwoCell', 'Suspension']:
        chosen_phi = [4.0, 6.0] if system == 'TwoCell' else [3.8991, 5.9986]
        for phi in chosen_phi:
            y += data[system][phi][1][i] # t_rot
            x += data[system][phi][0][i] # Cas
    model = np.polyfit(x, y, 1)
    #print('{}System'.format(system), 'doublet fitting line' if i else 'singlet fitting line', ':', model[0])
    print('doublet fitting line' if i else 'singlet fitting line', ':', model[0])
    predict = np.poly1d(model)
    x_fit = np.linspace(min(x), max(x), 100)
    y_fit = predict(x_fit)
    ax.plot(x_fit, y_fit, linestyle='--', linewidth=3, label='doublet fitting line' if i else 'singlet fitting line', color = 'r' if i else 'b')

ax.set_xlabel("Ca", fontsize = 20)
ax.set_ylabel("{}".format(r'$\gamma _{rot}$'), fontsize = 20)
ax.tick_params(labelsize = 15)
ax.xaxis.set_major_locator(MultipleLocator(0.05))
ax.xaxis.set_minor_locator(MultipleLocator(0.025))
ax.yaxis.set_major_locator(MultipleLocator(5))
ax.yaxis.set_minor_locator(MultipleLocator(2.5))
#ax.legend(frameon=False, bbox_to_anchor=(1.0, 1.0), loc='upper left')
#ax.legend(frameon=False)
#fig.tight_layout()

if SHOW:
    plt.show()
else:
    fig.savefig("Pictures/RotationTime_vs_Ca_combined_withoutLegend.png", dpi = 200)
    plt.close()

if SAVE:
    with open("T_rot_TCandS_w_{}_k_{}.pickle".format(window_width, threshold), 'wb') as handle:
        pickle.dump(data, handle, protocol=pickle.HIGHEST_PROTOCOL)

    print('Total time elapsed = {}'.format(str(datetime.timedelta(seconds=time.time()-start_time))))
# ===============================================================================
# Copyright 2021 An-Jun Liu
# Last Modified Date: 07/02/2021
# ===============================================================================
import numpy as np 
import ctypes
import logging
import sys
from scipy import stats
import math
import matplotlib.pyplot as plt 
import time
import datetime
import random

"""
This code is to make the system's roation time (seperate doublet and siglets state) vs Ca 
with varied volume fractions for suspension and two-cell system.
"""

path_preprocess = "/userdata4/ajliu/Data_Transfer/"
path_CT = "/raid6/ctliao/Data/HI_ordering/"
path_AJ = "/userdata4/ajliu/RBC_doublet/Data/"

def getRotationTime(input_phi, input_Ca, input_criteria_T, input_criteria_Dms, depend, system):

    # C extension
    lib = ctypes.cdll.LoadLibrary('./RBC_Utilities_CExtension.so')

    # Parameter
    # ===============================================================================
    phi = input_phi
    Ca = input_Ca
    Dm = 15.64
    #criteria_Ts = input_criteria_Ts if type(input_criteria_Ts) is list else [input_criteria_Ts]
    criteria_Dms = input_criteria_Dms if type(input_criteria_Dms) is list else [input_criteria_Dms]
    
    # Two-cell system
    if system == 0:
        angle = depend
        ncycle = 4000
        job_name = "phi_{}_Re_0.1_Ca_{}_aggregation_1KT_ncycle_{}_np_2_angle_{}".format(phi, Ca, ncycle, angle)
    
    # Suspension system
    else:
        ensemble_id = depend
        job_name = "h24phi{}Re0.1Ca{}WCA1zero0.8-{}".format(phi, Ca, ensemble_id)

    
    # Read the preprocessed files
    # ===============================================================================
    try:
        # Read the parameters
        with open(path_preprocess + "{}_parameter.txt".format(job_name)) as f:
            pre_parameters = f.readlines()
        timesteps = int((pre_parameters[pre_parameters.index("timesteps\n")+1])[:-1]) # number of COM files
        particle_numbers = int((pre_parameters[pre_parameters.index("particle_numbers\n")+1])[:-1])
        points_per_particle = int((pre_parameters[pre_parameters.index("points_per_particle\n")+1])[:-1])
        interval = int((pre_parameters[pre_parameters.index("interval\n")+1])[:-1]) # number of Ypos_t files
        
        dim = np.zeros(3, dtype = np.int32)
        tmp_dim = eval((pre_parameters[pre_parameters.index("dim\n")+1])[:-1])
        for i in range(3):
            dim[i] = tmp_dim[i]

        # Read the position of center of mass
        COMs = np.load(path_preprocess + "{}_COMs.npy".format(job_name))
        COMs_NB = np.load(path_preprocess + "{}_COMs_NB.npy".format(job_name))

        # Read the y positions
        Ypos_t = np.load(path_preprocess + "{}_Ypos_t.npy".format(job_name))
    
    except OSError:
        logging.error(" calcDoubletFraction():\nNo preprocessed data for simulation: {}\n".format(job_name))
        sys.exit(1)


    # Note: All the following computations can be parallelized

    # Calculate t_rot
    # ===============================================================================
    # Format of Ypos_t is Ypos_t[t, node_id]
    Periods = np.zeros(particle_numbers*points_per_particle)
    cutoff_frequency = 0.0135
    cutoff = int(cutoff_frequency*timesteps + 1) # igonre frequency lower than cutoff_frequency, i.e. rotation time > 1/cutoff_frequency starins

    for i in range(particle_numbers*points_per_particle):
        P = np.fft.rfft((Ypos_t[:, i] - np.mean(Ypos_t[:, i]))) # remove the DC term
        Periods[i] = interval/(np.argmax(np.abs(P[cutoff:]))+cutoff)

    # convert from the time unit in nodeposition format to COM format
    # in old output format, time unit in nodeposition format = WriteConfig and time unit in COM format = WriteProps
    # in new output format, time unit in nodeposition format = time unit in COM format = WriteProps
    # Updated Date: 01/29/2021 by An-Jun Liu
    scale = timesteps/interval
    rotation_time = stats.trim_mean(Periods, 0.1)*scale # remove the possible outliers



    # Calculate pair distance
    # ===============================================================================
    # Format of COMs is COMs[particle_id, t, 3]
    
    number_of_pairs = int((particle_numbers-1)*particle_numbers/2)
    diffpos = np.zeros((number_of_pairs, timesteps), dtype = np.float64)
    uncorrected_diffpos = np.zeros((number_of_pairs, timesteps), dtype = np.float64)
    indice_pairs = np.zeros((number_of_pairs, 2), dtype = np.int32)
    
    count = 0
    for i in range(particle_numbers-1):
        for j in range(i+1, particle_numbers):
            indice_pairs[count, 0], indice_pairs[count, 1] = i, j
            uncorrected_diffpos[count, :] = np.linalg.norm((COMs[i, :, :] - COMs[j, :, :]), axis=1)
            count += 1

    c_correctDiffpos = lib.correctDiffpos
    c_correctDiffpos(ctypes.c_void_p(diffpos.ctypes.data), ctypes.c_void_p(COMs.ctypes.data), ctypes.c_int(number_of_pairs),
    ctypes.c_int(timesteps), ctypes.c_double(Dm), ctypes.c_void_p(indice_pairs.ctypes.data), ctypes.c_void_p(dim.ctypes.data)) 
    

          
    # Calculate doublet fraction
    # ===============================================================================
    criteria_T = input_criteria_T
    period = int(round(criteria_T*rotation_time))
    DF_indivisual = np.zeros((particle_numbers, timesteps-period))
    window_width = 100
    threshold = 0.9
    target1, target2 = int(window_width*threshold), int(window_width*(1-threshold))

    for criteria_Dm in criteria_Dms:
        doublet_or_not = np.zeros((number_of_pairs, timesteps - period), dtype = np.int32) # 1 means there is a doublet
        state_series = np.zeros((number_of_pairs, timesteps - period), dtype = np.int32)
        end_time = np.zeros(1, dtype = np.int32)

        c_calcDF = lib.calcDF

        c_calcDF(ctypes.c_void_p(doublet_or_not.ctypes.data), ctypes.c_void_p(diffpos.ctypes.data), ctypes.c_int(period),
        ctypes.c_int(timesteps), ctypes.c_int(number_of_pairs), ctypes.c_double(Dm), ctypes.c_double(criteria_Dm),
        ctypes.c_void_p(COMs_NB.ctypes.data), ctypes.c_void_p(indice_pairs.ctypes.data), ctypes.c_void_p(state_series.ctypes.data),
        ctypes.c_void_p(end_time.ctypes.data), ctypes.c_void_p(dim.ctypes.data))

        if np.average(doublet_or_not) < 0.5:
            return False
        else:
            print("phi = {}, Ca = {}, angle = {}".format(input_phi, input_Ca, depend))

        # recognize the state time series for each particle
        for i in range(number_of_pairs):
            j, k = indice_pairs[i, 0], indice_pairs[i, 1]

            # go through the i-th pair
            for t in range(timesteps-period):
                if doublet_or_not[i, t] > 0:
                    DF_indivisual[j, t] = 1
                    DF_indivisual[k, t] = 1

        # calculate the rotation time of doublet and singlet state for each particle
        t_rot_stat = [[], []] # [[singlet], [doublet]]

        for i in range(particle_numbers):
            # split the time series according to the state
            # initialize
            split = []
            record = 0
            current_state = 0

            df = np.sum(DF_indivisual[i, :window_width])
            if df > target1:
                record = 1
                current_state = 1
                st = 0
            elif df < target2:
                record = 1
                current_state = 0
                st = 0

            # go throught the tim series:
            for t in range(window_width, timesteps-period-window_width):
                df += (DF_indivisual[i, t] - DF_indivisual[i, t-window_width]) # update DF
                if record:
                    if (current_state == 1 and df < target1) or (current_state == 0 and df > target2):
                        if (t - st) >= window_width:
                            split.append([st, t, current_state])
                        record = 0
                else:
                    if df >= target1:
                        st = t
                        current_state = 1
                        record = 1
                    elif df <= target2:
                        st = t
                        current_state = 0
                        record = 1

            # extreme case
            if record and (t - st) >= window_width:
                split.append([st, t, current_state])


            # calculate rotation time
            for slice in split:
                for j in range(points_per_particle):
                    cutoff = int(cutoff_frequency*(slice[1]-slice[0]) + 1)
                    P = np.fft.rfft((Ypos_t[int(slice[0]/scale):int(slice[1]/scale), i] - np.mean(Ypos_t[int(slice[0]/scale):int(slice[1]/scale), i]))) # remove the DC term
                    t_rot_stat[slice[2]].append((slice[1]-slice[0])/(np.argmax(np.abs(P[cutoff:]))+cutoff))

                if slice[2]:
                    fig, (ax1, ax2) = plt.subplots(1, 2, figsize = (16, 7))
                    st, et = slice[0], slice[1]
                    for j in range(points_per_particle):
                        ax1.plot(np.arange(int(slice[0]/scale), int(slice[1]/scale)), 
                        (Ypos_t[int(slice[0]/scale):int(slice[1]/scale), i] - np.mean(Ypos_t[int(slice[0]/scale):int(slice[1]/scale), i]))*0.5, label = "node {}".format(j))
                        ax1.set_xlabel("time ({})".format(r'$\dot \gamma t$'), fontsize = 20)
                        ax1.set_ylabel("y ({})".format(r'$\mu m$'), fontsize = 20)
                        ax1.legend(prop={'size': 12})
                        ax1.set_title("y(t)", fontsize = 20)
                        ax1.xaxis.set_tick_params(labelsize = 12)
                        ax1.yaxis.set_tick_params(labelsize = 12)


                    for j in range(points_per_particle):
                        #yf = np.fft.rfft(data[:, indices[i]])
                        yf = np.fft.rfft(data[:, i])
                        xf = np.fft.rfftfreq(N, d = 1.0/fs)
                        ax2.plot(xf[offset:], (2.0/N * np.abs(yf))[offset:], label = "node {}".format(i+1))
                    ax2.set_xlabel("frequency ({})".format(r'$\dfrac{1}{\dot \gamma t}$'), fontsize = 20)
                    ax2.set_ylabel("Amplitude", fontsize = 20)
                    ax2.xaxis.set_tick_params(labelsize = 12)
                    ax2.yaxis.set_tick_params(labelsize = 12)
                    ax2.legend(prop={'size': 12})
                    plt.title("Corresponding Frequency Spectrum", fontsize = 20)
                    fig.tight_layout()
                    plt.show()
        
        return True



### main test code
start_time = time.time()
### Two Cell

vol = 746.3163
phi_range =[]
for x in range(30, 41):
    phi = 2*vol/(24*x**2)
    phi_range.append(float(str(phi*100)[:3]))
Ca_range = [i*0.01 for i in range(1, 21)]
angle_range = [90 - 10*j for j in range(18)]
phi_range.sort()

# find valid set of parameters
st = time.time()
found = False
while not found and (time.time() - st) < 5*60:
    phi = random.choice(phi_range)
    Ca = random.choice(Ca_range[14:])
    angle = random.choice(angle_range)
    found = getRotationTime(phi, Ca, 1, 1, angle, 0)

print('Total time elapsed = {}'.format(str(datetime.timedelta(seconds=time.time()-start_time))))
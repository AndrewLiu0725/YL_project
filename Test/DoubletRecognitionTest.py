# ===============================================================================
# Copyright 2021 An-Jun Liu
# Last Modified Date: 06/14/2021
# ===============================================================================
import numpy as np
import matplotlib.pyplot as plt
import sys
import time
from scipy import stats
import ctypes
import os
import logging
import traceback

"""
This utility python file provides functions to calculate doublet fraction, interparticle/elastic stress, intrinsic viscosity, 
and relative viscosity for two-cell and suspension system. 
The Multi-Line Docstring for each function is also provided. 
Please refer to these docstrings before using these utility functions.
"""


path_preprocess = "data/"


def calcDoubletFraction(input_phi, input_Ca, input_criteria_T, input_criteria_Dms, depend, system, st, et):
    """
    Input:
    
    input_criteria_Dms must be either lists or floats
    
    make_plot = 0 means doesn't need to make plot, 1 means saving df vs t plot without playing it, 2 means displaying the plot

    system = 0 means two-cell system, 1 means suspension

    depend means angle when system = 0, esemble id when system = 1
    
    Output:
    
    [[doublet fraction], timesteps, [state series], diffpos, particle_numbers]
    
    In state series, 1 means doublet, 2 means kayaking, 3 means sliding, and 4 means transition states

    Exit code:
    0 Success
    1 OSError

    Note:
    1. All the arrays shared by python and C program are stored in double so that we don't need to consider the data compatibility.
    2. You have to preprocess the data before using this utility function.
    """

    # C extension
    lib = ctypes.cdll.LoadLibrary('./RBC_Utilities_CExtension.so')

    # Parameter
    # ===============================================================================
    phi = input_phi
    Ca = input_Ca
    Dm = 15.64
    criteria_Dm = input_criteria_Dms
    
    # Two-cell system
    if system == 0:
        angle = depend
        ncycle = 4000
        job_name = "phi_{}_Re_0.1_Ca_{}_aggregation_1KT_ncycle_{}_np_2_angle_{}".format(phi, Ca, ncycle, angle)
    
    # Suspension system
    elif system == 1:
        ensemble_id = depend
        job_name = "h24phi{}Re0.1Ca{}WCA1zero0.8-{}".format(phi, Ca, ensemble_id)

    else:
        logging.error(" calcDoubletFraction():\nWrong system code! The acceptable value is 0 or 1.")
        sys.exit(1)

    
    # Read the preprocessed files
    # ===============================================================================
    try:
        # Read the parameters
        with open(path_preprocess + "{}_parameter.txt".format(job_name)) as f:
            pre_parameters = f.readlines()
        timesteps = int((pre_parameters[pre_parameters.index("timesteps\n")+1])[:-1]) # number of COM files
        particle_numbers = int((pre_parameters[pre_parameters.index("particle_numbers\n")+1])[:-1])
        points_per_particle = int((pre_parameters[pre_parameters.index("points_per_particle\n")+1])[:-1])
        #interval = int((pre_parameters[pre_parameters.index("interval\n")+1])[:-1]) # number of Ypos_t files
        
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
        P = np.fft.rfft((Ypos_t[:, i] - np.mean(Ypos_t[:, i]))) # remove the DC term to eliminate the large amplitude 0 Hz component
        Periods[i] = timesteps/(np.argmax(np.abs(P[cutoff:]))+cutoff)
        '''
        math:
        f = [0, 1/nd, ..., (n/2)/nd] where nd is fixed as the real time interval 
        and f's unit is also in real time (i.e. strains).
        Thus, f_cutoff = i_cutoff / nd -> i_cutoff = f_cutoff x timesteps
        Moreover, period = 1/f = 1/((i_max+i_cutoff)/nd) = timesteps / (i_max+i_cutoff)
        '''


    rotation_time = stats.trim_mean(Periods, 0.1) # remove the possible outliers
    print("rotation time =", rotation_time)



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
    fig, ax1 = plt.subplots(figsize = (8, 6))

    output_DF, output_state = [], []
    criteria_T = input_criteria_T
    period = int(round(criteria_T*rotation_time))

    doublet_or_not = np.zeros((number_of_pairs, timesteps - period), dtype = np.int32) # 1 means there is a doublet
    state_series = np.zeros((number_of_pairs, timesteps - period), dtype = np.int32)
    end_time = np.zeros(1, dtype = np.int32)

    c_calcDF = lib.calcDF

    c_calcDF(ctypes.c_void_p(doublet_or_not.ctypes.data), ctypes.c_void_p(diffpos.ctypes.data), ctypes.c_void_p(uncorrected_diffpos.ctypes.data),
    ctypes.c_void_p(COMs.ctypes.data), ctypes.c_void_p(COMs_NB.ctypes.data), ctypes.c_void_p(indice_pairs.ctypes.data), ctypes.c_void_p(dim.ctypes.data),
    ctypes.c_int(period), ctypes.c_int(timesteps), ctypes.c_int(number_of_pairs), ctypes.c_double(Dm), ctypes.c_double(criteria_Dm), ctypes.c_void_p(end_time.ctypes.data))
    
    # Assume that there is no aggregation structure which consists more than 2 partilces
    number_of_doublets = np.sum(doublet_or_not, axis = 0)
    # one doublet has two RBCs so to calculate doublet fraction, we need multiply # of doublets by 2
    output_DF.append(number_of_doublets*2/particle_numbers)
    output_state.append(state_series)
    
    ax1.plot(np.arange(st, et), number_of_doublets[st: et], color = 'b', label = "DF (min Max algo)")
        


    # Make plot
    # ===============================================================================
    ax1.set_xlabel("timesteps (in WriteProps unit)", fontsize = 20)
    ax1.set_ylabel("doublet fraction", fontsize = 20)
    if system == 0:
        ax1.set_title("Two-cell system\nphi={}, Ca={}, angle={}, Re=0.1, h=24".format(phi, Ca, angle), fontsize = 20)
    else:
        ax1.set_title("Suspension\nphi={}, Ca={}, esemble {}\nRe=0.1, D=1, eqWCA=0.8, h=24".format(phi, Ca, ensemble_id), fontsize = 20)
    ax2 = ax1.twinx()
    ax2.set_ylabel('distance', fontsize = 20)
    ax2.plot(np.arange(st, et), diffpos[0, st:et], color = 'r', label = "corrected diff")
    ax2.plot(np.arange(st, et), uncorrected_diffpos[0, st:et], color = 'g', label = "uncorrected diff")
    h1, l1 = ax1.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax1.legend(h1+h2, l1+l2, loc=0)
    fig.tight_layout()
    plt.show()

    #return [output_DF, timesteps, output_state, diffpos, COMs, COMs_NB]
    return [output_DF, end_time[0], output_state, diffpos, particle_numbers]


# suspension system
#[2.3994, 2.9993, 3.4492, 3.8991, 4.9488, 5.9986]
#calcDoubletFraction(3.4492, 0.01, 1, 1, 2, 5, 1)

# two-cell system
#calcDoubletFraction(6.0, 0.06, 1, 1, -80, 0, 1700, 1900)
#calcDoubletFraction(3.8, 0.02, 1, 1, 40, 0, 0, 400)
#calcDoubletFraction(4.7, 0.05, 1, 1, 0, 0, 0, 400)
calcDoubletFraction(4.7, 0.15, 1, 1, 30, 0, 0, 500)
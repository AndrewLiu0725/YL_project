# ===============================================================================
# Copyright 2021 An-Jun Liu
# Last Modified Date: 06/10/2021
# ===============================================================================
import numpy as np 
import ctypes
import logging
import sys
from scipy import stats
import matplotlib.pyplot as plt 
from RBC_Utilities import getSuspensionParameterSets

"""
This code is to make the system's roation time vs Ca with varied volume fractions for suspension and two-cell system.
"""

path_preprocess = "/userdata4/ajliu/Data_Transfer/"
path_CT = "/raid6/ctliao/Data/HI_ordering/"
path_AJ = "/userdata4/ajliu/RBC_doublet/Data/"

def getRotationRate(input_phi, input_Ca, input_criteria_T, input_criteria_Dms, make_plot, depend, system):
    # C extension
    lib = ctypes.cdll.LoadLibrary('./cforDoublet_Functions.so')

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

    offset = int(0.0125*timesteps)
    for i in range(particle_numbers*points_per_particle):
        P = np.fft.rfft(Ypos_t[:, i])
        Periods[i] = interval/(np.argmax(np.abs(P[offset:]))+offset) # divide timesteps by 2 since here is the number of the timesteps of bond0.vtk

    # convert from the time unit in nodeposition format to COM format
    # in old output format, time unit in nodeposition format = WriteConfig and time unit in COM format = WriteProps
    # in new output format, time unit in nodeposition format = time unit in COM format = WriteProps
    # Updated Date: 01/29/2021 by An-Jun Liu
    rotation_time = stats.trim_mean(Periods, 0.1)*(timesteps/interval) 
    #rotation_time = np.mean(Periods)*(timesteps/interval) 

    return rotation_time


### main code
### Two Cell
vol = 746.3163
phi_range =[]
for x in range(30, 41):
    phi = 2*vol/(24*x**2)
    phi_range.append(float(str(phi*100)[:3]))
Ca_range = [i*0.01 for i in range(1, 21)]
angle_range = [90 - 10*j for j in range(18)]


fig, ax = plt.subplots(figsize = (8, 6))

for phi_index, phi in enumerate(phi_range):
    if phi_index % 2 != 0: continue
    t_rot = np.zeros(len(Ca_range))
    for Ca_index, Ca in enumerate(Ca_range):
        for angle_index, angle in enumerate(angle_range):
            t_rot[Ca_index] += getRotationRate(phi, Ca, 1, 1, 0, angle, 0)
    ax.plot(Ca_range, t_rot/len(angle_range), label = '{} = {}%'.format(r'$\phi$', round(phi,1)))

ax.set_title("Rotation Time vs Ca (Two-Cell System)", fontsize = 20)
ax.set_xlabel("Ca", fontsize = 12)
ax.set_ylabel("Rotation Time ({})".format(r'$\dot \gamma t$'), fontsize = 12)
ax.legend(fontsize = 12)
plt.savefig("Pictures/TwoCellSystem_RotationTime_vs_Ca.png", dpi = 200)
plt.close()

### Suspension
fig, ax = plt.subplots(figsize = (8, 6))
[phis, parameter_set] = getSuspensionParameterSets()
phis.sort()
for phi in phis:
    t_rot = np.zeros(len(list(parameter_set[phi].keys())))
    Ca_range = list(parameter_set[phi].keys())
    Ca_range.sort()
    for Ca_index, Ca in enumerate(Ca_range):
        ensemble_count = 0
        for ensemble_id in parameter_set[phi][Ca]:
            try:
                t_rot[Ca_index] += getRotationRate(phi, Ca, 1, 1, 0, ensemble_id, 1)
                ensemble_count += 1
            except:
                pass
        t_rot[Ca_index] = t_rot[Ca_index]/ensemble_count
    ax.plot(Ca_range, t_rot*4000/3669, label = '{} = {}%'.format(r'$\phi$', round(phi,1)))

ax.set_title("Rotation Time vs Ca (Suspension System)", fontsize = 20)
ax.set_xlabel("Ca", fontsize = 12)
ax.set_ylabel("Rotation Time ({})".format(r'$\dot \gamma t$'), fontsize = 12)
ax.legend(fontsize = 12)
plt.savefig("Pictures/SuspensionSystem_RotationTime_vs_Ca.png", dpi = 200)
plt.close()
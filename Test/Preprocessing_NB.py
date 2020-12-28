import os
import sys
import numpy as np 
import time

total_start_time = time.time()
start_time = time.time()

# parameter to be assigned, input format as mode_{}_phi_{}_Ca_{}_D_{}_eqWCA_{} / mode_{}_phi_{}_Ca_{}_angle_{}
###########################################################################
parameters = sys.argv[1].split('_')
phi = parameters[parameters.index('phi')+1]
Ca = parameters[parameters.index('Ca')+1]
mode = parameters[parameters.index('mode')+1]

if mode == '0':
    # Suspension system
    particle_num_dict = {'5': 33, '4': 26}
    particle_numbers = particle_num_dict[phi]

    D = parameters[parameters.index('D')+1]
    eqWCA = parameters[parameters.index('eqWCA')+1]
    if phi == "5":
        job_name = "phi{}/Re0.1/eqWCA{}/h24phi{}Re0.1Ca{}D{}eqWCA{}".format(phi, eqWCA, phi, Ca, D, eqWCA)
    else:
        job_name = "phi{}/eqWCA{}/h24phi{}Re0.1Ca{}D{}eqWCA{}".format(phi, eqWCA, phi, Ca, D, eqWCA)
    path_job = "/userdata4/ctliao/Project/HI_ordering/{}/data".format(job_name)
    dim = [144, 24, 144]
    filename_prefix = "/userdata4/ajliu/Data_Transfer/h24phi{}Re0.1Ca{}D{}eqWCA{}".format(phi, Ca, D, eqWCA)

    f = open(path_job+'/sphere_props.0.dat')
    data = f.readlines()
    timesteps = len(data) - 1
    f.close()


else:
    # Two-cell system 
    particle_numbers = 2

    angle = parameters[parameters.index('angle')+1]
    job_name = "phi_{}_Re_0.1_Ca_{}_aggregation_1KT_ncycle_2000_np_2_angle_{}".format(phi, Ca, angle)
    path_job = "/userdata4/ajliu/RBC_doublet/Data/{}/data".format(job_name)
    path_parameter = "/userdata4/ajliu/RBC_doublet/parameter/{}.dat".format(job_name)

    fp = open(path_parameter)
    dim = [int(axis) for axis in fp.readlines()[3].split()]
    fp.close()

    filename_prefix = "/userdata4/ajliu/Data_Transfer/{}".format(job_name)

    f = open(path_job+'/sphere_props.0.dat')
    data = f.readlines()
    timesteps = len(data) - 1
    f.close()



# Get COM data here
###########################################################################
def getCOM(path):
    # shape of temp_positionCOM would be [timesteps][3]
    COM = np.zeros((timesteps, 3))
    f = open(path)
    data = f.readlines()
    for t in range(timesteps):
        COM[t, :] = [float(i) for i in data[t+1].split()[1: 4]]
    f.close()
    return COM

COMs = np.zeros((particle_numbers, timesteps, 3), dtype = np.float64) # stored in double

for i in range(particle_numbers):
    COMs[i, :, :] = getCOM(path_job+'/sphere_props.{}.dat'.format(i))

# Write COMs, bond0, and parameter
###########################################################################
np.save(filename_prefix + "_COMs_NB.npy", COMs)
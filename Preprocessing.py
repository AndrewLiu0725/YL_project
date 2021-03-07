# ===============================================================================
# Copyright 2021 An-Jun Liu
# Last Modified Date: 03/07/2021
# ===============================================================================
import os
import sys
import numpy as np 
import time
import logging

"""
This file is to preprocess the data for calculating the doublet function.
The preprocessed data includes
1. time series of bounded center of mass of each particle, i.e. _COMs.npy
2. time series of unbounded center of mass of each particle, i.e. _COMs_NB.npy
3. time series of y coordinate of specific nodes, i.e. _Ypos_t.npy
4. parameters including timesteps (COMs time unit), interval (Ypos time unit), particle numbers,
points per particle (used in Ypos_t), and dimensions, i.e. _parameter.txt

Usage:
system=0 means two-cell system, 1 means suspension.
verbose=0 means only showing error meassage, 1 means also showing warnning message
e.g. python3 Preprocessing.py [simulation folder name] system=0 verbose=0

Exit code:
0 success
1 already run
2 OSError
3 IndexError
4 Wrong command format
"""

def BinarySearch(prefix, expected_end, time_increment):
    """
    Example of prefix:
    
    "/raid6/ctliao/Data/HI_ordering/h24_phi3.8991_Re0.1_Ca0.07_WCA1_zero0.8-7/data/nodePositions_.dat"

    Output:
    
    number of valid files
    """

    if (os.path.isfile(prefix.format((expected_end-1)*time_increment))): # complete simulation
        return expected_end

    elif not (os.path.isfile(prefix.format(0))): # no base file
        logging.error(" No such file: {}".format(prefix.format(0)))
        sys.exit(2)

    else: # binary search here
        L = 0
        R = expected_end - 1
        mid = int((L+R)/2)
        while(L != mid):
            if (os.path.isfile(prefix.format(mid*time_increment))):
                L = mid
            else:
                R = mid
            mid = int((L+R)/2)
        return mid+1

if (len(sys.argv) < 4):
    logging.error(" Wrong command format!\nUsage:\npython3 Preprocessing.py [simulation folder name] system=[0 or 1] verbose=[0 or 1]")
    sys.exit(4)


# Setup
# ===============================================================================
CT_folder = "/raid6/ctliao/Data/HI_ordering/"
AJ_folder = "/userdata4/ajliu/RBC_doublet/Data/"
bead_number = 642 # 642 beads per particle
points_per_particle = 6
job_name = sys.argv[1] # name of the data folder

# catch flags
system_found, verbose_found = False, False
for i in range(2, len(sys.argv)):
    if len(sys.argv[i]) == 8 and sys.argv[i][:6] == "system":
        system = int(sys.argv[i][-1])
        system_found = True
    elif len(sys.argv[i]) == 9 and sys.argv[i][:7] == "verbose":
        verbose = int(sys.argv[i][-1])
        verbose_found = True
if (not system_found) or (not verbose_found):
    logging.error(" Wrong command format!\nUsage:\npython3 Preprocessing.py [simulation folder name] system=[0 or 1] verbose=[0 or 1]")
    sys.exit(4)


# Two-cell system 
# job_name format is phi_5.0_Re_0.1_Ca_0.2_aggregation_1KT_ncycle_2000_np_2_angle_70
if system == 0:
    working_folder = AJ_folder
    tmp_filename_prefix = job_name

# Suspension system
# job_name format is either h24_phi4.4989_Re0.1_Ca0.06_WCA1_zero0.8-8 or h24phi4.9488Re0.1Ca0.06WCA1zero0.8-4
# tmp_filename_prefix format would be like this: h24phi4.9488Re0.1Ca0.06WCA1zero0.8-4
# Updated Date: 01/29/2021 by An-Jun Liu
else:
    working_folder = CT_folder
    tmp_filename_prefix = ''.join(job_name.split('_'))

path_job = working_folder + job_name + "/data/"
filename_prefix = "/userdata4/ajliu/Data_Transfer/" + tmp_filename_prefix # filename prefix of the preprocessed data

# check data format
if os.path.isfile(path_job+"sphere_props.0.dat"):
    new_data_format_flag = 0 # 0 means old data format, 1 means new one
else:
    if os.path.isfile(path_job+"nodePositions0.dat"):
        new_data_format_flag = 1
    else:
        logging.error(" Neither {} nor {} exists!".format(path_job+"sphere_props.0.dat", path_job+"nodePositions0.dat"))
        sys.exit(2)



# Get the parameters (may need more robustic way to read the parameters)
# ===============================================================================
try:
    f = open(working_folder + job_name + "/init/parameter.dat", 'r')
    data = f.readlines()
    dim = [int(axis) for axis in data[3].split()]
    particle_numbers = int(data[1].split()[0])
    WriteProps = int(data[11].split()[0])
    WriteConfig = int(data[11].split()[1])
    expected_timesteps = int(int(data[5].split()[0])*int(data[5].split()[1])/WriteProps) # NumCycle x NumStep / WriteProps
    expected_interval = int(int(data[5].split()[0])*int(data[5].split()[1])/WriteConfig) # NumCycle x NumStep / WriteConfig
    f.close()

except OSError:
    logging.error(" No such file: {}".format(working_folder + job_name + "/init/parameter.dat"))
    sys.exit(2)

except IndexError:
    logging.error(" list index out of range in reading parameters from {}.\nMay need change the format of parameter.dat or the way this script reading parameters".format(working_folder + job_name + "/init/parameter.dat"))
    sys.exit(3)



# Check if this case is already preprocessed
# ===============================================================================
if os.path.isfile(filename_prefix + "_parameter.txt"):
    f = open(filename_prefix + "_parameter.txt", 'r')
    pre_parameters = f.readlines()
    previous_timesteps = int((pre_parameters[pre_parameters.index("timesteps\n")+1])[:-1])
    f.close()

    if new_data_format_flag:
        current_timesteps = BinarySearch(path_job+"nodePositions{}.dat", expected_timesteps, WriteProps)
    else:
        current_timesteps = sum(1 for line in open(path_job+'sphere_props.0.dat'))-1
    
    if current_timesteps <= previous_timesteps:
        if verbose:
            logging.warning(" Job \"{}\" is already preprocessed.".format(job_name))
        sys.exit(1)



# Old data format
# has sphere_props.{}.dat written in frequency WriteProps
# has bond0_t{}.vtk written in frequency WriteConfig
# ===============================================================================
def getCOM(path):
    # output shape would be [timesteps][3]
    COM = np.zeros((2, timesteps, 3))
    f = open(path, 'r')
    data = f.readlines()
    for t in range(timesteps):
        COM[0, t, :] = [(float(i)%dim[dim_index]) for dim_index, i in enumerate(data[t+1].split()[1: 4])]
        COM[1, t, :] = [float(i) for i in data[t+1].split()[1: 4]]
    f.close()
    return COM

point_offset = int(bead_number/points_per_particle)
def getYpos(time):
    f = open(path_job+"bond0_t{}.vtk".format(int(time) if time != 0 else "0000"), 'r')
    data = f.readlines()
    Ypos = []
    for i in range(particle_numbers):
        for j in range(points_per_particle):
            Ypos.append(float(data[i*bead_number + j*point_offset + 6].split()[1])) # 1 means y pos
    f.close()
    return np.array(Ypos)

if new_data_format_flag == 0:
    # Get COM data here
    # ===============================================================================
    timesteps = sum(1 for line in open(path_job+'sphere_props.0.dat'))-1
    COMs = np.zeros((particle_numbers, timesteps, 3), dtype = np.float64) # stored in double
    COMs_NB = np.zeros((particle_numbers, timesteps, 3), dtype = np.float64) # stored in double
    for i in range(particle_numbers):
        result = getCOM(path_job+'sphere_props.{}.dat'.format(i))
        COMs[i, :, :] = result[0]
        COMs_NB[i, :, :] = result[1]

    # Get bond0.vtk data here
    # ===============================================================================
    interval = BinarySearch(path_job+"bond0_t{}.vtk", expected_interval, WriteConfig)
    Ypos_t = np.zeros((interval, particle_numbers*points_per_particle))
    for i in range(interval):
        Ypos_t[i, :] = getYpos(int(i*WriteConfig))


# New data format 
# has nodePositions{}.dat written in frequency WriteProps
# has bond0_t{}.vtk written in frequency WriteConfig
# ===============================================================================
else:
    timesteps = BinarySearch(path_job+"nodePositions{}.dat", expected_timesteps, WriteProps)
    interval = timesteps
    increment = int(bead_number/points_per_particle)
    Ypos_t = np.zeros((timesteps, particle_numbers*points_per_particle))
    COMs = np.zeros((particle_numbers, timesteps, 3), dtype = np.float64) # stored in double
    COMs_NB = np.zeros((particle_numbers, timesteps, 3), dtype = np.float64) # stored in double
    for i in range(timesteps):
        t = i*WriteProps
        data = np.loadtxt(path_job+"nodePositions{}.dat".format(t))
        # Ypos_t
        for j in range(particle_numbers*points_per_particle):
            Ypos_t[i, j] = data[j*increment, 1]

        # COMS_NB and COMs_NB
        for j in range(particle_numbers):
            COMs_NB[j, i, :] = np.sum(data[j*bead_number:(j+1)*bead_number, :], axis = 0)/bead_number
            COMs[j, i, 1] = COMs_NB[j, i, 1]
            COMs[j, i, 0] = np.sum(np.mod(data[j*bead_number:(j+1)*bead_number, 0], dim[0]))/bead_number
            COMs[j, i, 2] = np.sum(np.mod(data[j*bead_number:(j+1)*bead_number, 2], dim[2]))/bead_number

if (timesteps < expected_timesteps) and verbose:
    logging.warning(" Job \"{}\" is incomplete!".format(job_name))



# Write COMs, bond0, and parameter into files
# ===============================================================================
np.save(filename_prefix + "_COMs.npy", COMs)
np.save(filename_prefix + "_COMs_NB.npy", COMs_NB)
np.save(filename_prefix + "_Ypos_t.npy", Ypos_t)
f = open(filename_prefix+"_parameter.txt", "w")
f.write("timesteps\n{}\n".format(timesteps))
f.write("particle_numbers\n{}\n".format(particle_numbers))
f.write("interval\n{}\n".format(interval))
f.write("points_per_particle\n{}\n".format(points_per_particle))
f.write("dim\n{}\n".format(dim))
#f.write("WriteProps\n{}\nWriteConfig\n{}\n".format(WriteProps, WriteConfig))
f.close()
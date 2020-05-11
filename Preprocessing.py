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
    timesteps = 2500

    D = parameters[parameters.index('D')+1]
    eqWCA = parameters[parameters.index('eqWCA')+1]
    if phi == "5":
        job_name = job_name = "phi{}/Re0.1/eqWCA{}/h24phi{}Re0.1Ca{}D{}eqWCA{}".format(phi, eqWCA, phi, Ca, D, eqWCA)
    else:
        job_name = "phi{}/eqWCA{}/h24phi{}Re0.1Ca{}D{}eqWCA{}".format(phi, eqWCA, phi, Ca, D, eqWCA)
    path_job = "/userdata4/ctliao/Project/HI_ordering/{}/data".format(job_name)
    dim = [144, 24, 144]
    filename_prefix = "/userdata4/ajliu/Data_Transfer/h24phi{}Re0.1Ca{}D{}eqWCA{}".format(phi, Ca, D, eqWCA)


else:
    # Two-cell system 
    particle_numbers = 2
    timesteps = 2000

    angle = parameters[parameters.index('angle')+1]
    job_name = "phi_{}_Re_0.1_Ca_{}_aggregation_1KT_ncycle_2000_np_2_angle_{}".format(phi, Ca, angle)
    path_job = "/userdata4/ajliu/RBC_doublet/Data/{}/data".format(job_name)
    path_parameter = "/userdata4/ajliu/RBC_doublet/parameter/{}.dat".format(job_name)

    fp = open(path_parameter)
    dim = [int(axis) for axis in fp.readlines()[3].split()]
    fp.close()

    filename_prefix = "/userdata4/ajliu/Data_Transfer/{}".format(job_name)



# Get COM data here
###########################################################################
def getCOM(path):
    # shape of temp_positionCOM would be [timesteps][3]
    COM = np.zeros((timesteps, 3))
    f = open(path)
    data = f.readlines()
    for t in range(timesteps):
        COM[t, :] = [(float(i)%dim[dim_index]) for dim_index, i in enumerate(data[t+1].split()[1: 4])]
    return COM

COMs = np.zeros((particle_numbers, timesteps, 3), dtype = np.float64) # stored in double

for i in range(particle_numbers):
    COMs[i, :, :] = getCOM(path_job+'/sphere_props.{}.dat'.format(i))
print('Time elpased to collect COM data = ', time.time()-start_time)



# Get bond0.vtk data here
###########################################################################
start_time = time.time()

bead_number = 642 # 642 beads per particle
points_per_particle = 6
time_index = [] # every WriteConfig timesteps will generate a bond0.vtk, len(time_index) = timesteps/2
for fn in os.listdir(path_job):
    if fn.split('_')[0] == "bond0": time_index.append(((fn.split('_')[1]).split('.')[0])[1:])
sorted(time_index, key= lambda k: int(k))

def getYpos(time):
    f = open(path_job+"/bond0_t{}.vtk".format(time))
    data = f.readlines()
    entering_data = 0
    for line_index, line in enumerate(data):
        if line[0].isdigit() and entering_data == 0:
            end_pos_format = line_index
            entering_data = 1
            break
    Ypos = []
    for i in range(particle_numbers):
        for j in range(points_per_particle):
            Ypos.append(float(data[i*bead_number + j*int(bead_number/points_per_particle) + end_pos_format].split()[1])) # 1 means y pos
    f.close()
    return Ypos


interval = int(timesteps/2)
Ypos_t = []
for i in range(interval):
    Ypos_t.append(getYpos(time_index[i]))
Ypos_t = np.array(Ypos_t)
print('Time elpased to collect bond0.vtk data = ', time.time()-start_time)

# Write COMs, bond0, and parameter
###########################################################################
np.save(filename_prefix + "_COMs.npy", COMs)
np.save(filename_prefix + "_Ypos_t.npy", Ypos_t)

f = open(filename_prefix+"_parameter.txt", "w")
f.write("timesteps\n{}\n".format(timesteps))
f.write("particle_numbers\n{}\n".format(particle_numbers))
f.write("interval\n{}\n".format(interval))
f.write("points_per_particle\n{}\n".format(points_per_particle))
f.write("dim\n{}\n".format(dim))
f.close()
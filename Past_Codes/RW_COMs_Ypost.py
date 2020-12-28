import os
import sys
import time

total_start_time = time.time()
start_time = time.time()

# parameter to be assigned, input format as phi{}Re0.1Ca{}D{}eqWCA{}
parameters = sys.argv[1].split('_')
phi = parameters[parameters.index('phi')+1]
Ca = parameters[parameters.index('Ca')+1]
D = parameters[parameters.index('D')+1]
eqWCA = parameters[parameters.index('eqWCA')+1]

dim = [144, 24, 144]

# Get COM data here
def getCOM(path):
    # shape of temp_positionCOM would be [timesteps][3]
    temp_positionCOM = []
    with open(path) as f:
        for index, line in enumerate(f):
            if index > 0: temp_positionCOM.append([(float(i)%dim[dim_index]) for dim_index, i in enumerate(line.split()[1: 4])])
    return temp_positionCOM

if phi == "5":
    job_name = job_name = "phi{}/Re0.1/eqWCA{}/h24phi{}Re0.1Ca{}D{}eqWCA{}".format(phi, eqWCA, phi, Ca, D, eqWCA)
else:
    job_name = "phi{}/eqWCA{}/h24phi{}Re0.1Ca{}D{}eqWCA{}".format(phi, eqWCA, phi, Ca, D, eqWCA)
path_job = "/userdata4/ctliao/Project/HI_ordering/{}/data".format(job_name)

particle_numbers = 0
for fn in os.listdir(path_job):
    if fn.split('.')[0] == "sphere_props": particle_numbers += 1
print("number of particles =", particle_numbers)

timesteps = len(getCOM(path_job+'/sphere_props.0.dat')) # need to store

COMs = []
for i in range(particle_numbers):
    COMs.append(getCOM(path_job+'/sphere_props.{}.dat'.format(i)))
print('Time elpased to collect COM data = ', time.time()-start_time)


# Get bond0.vtk data here
start_time = time.time()

bead_number = 642 # 642 beads per particle
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
        Ypos.append(float(data[i*bead_number+end_pos_format].split()[1]))
    f.close()
    return(Ypos)


interval = 200
Ypos_t = []
for i in range(interval):
    Ypos_t.append(getYpos(time_index[int(timesteps/4)+i])) # from the middle of bond0.vtks
    
print('Time elpased to collect bond0.vtk data = ', time.time()-start_time)


# Write COMs, bond0, and parameter
filename_prefix = "h24phi{}Re0.1Ca{}D{}eqWCA{}".format(phi, Ca, D, eqWCA)

f = open(filename_prefix + "_COMs.txt", "w")
f.write(str(COMs))
f.close()

f = open(filename_prefix+"_Ypos_t.txt", "w")
f.write(str(Ypos_t))
f.close()

f = open(filename_prefix+"_parameter.txt", "w")
f.write("timesteps = {}\n".format(timesteps))
f.write("particle_numbers  ={}\n".format(particle_numbers))
f.close()
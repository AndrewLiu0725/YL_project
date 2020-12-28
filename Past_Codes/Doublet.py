import numpy as np 
import matplotlib.pyplot as plt
import os
import sys

# %(1*int(dim[dim_index]))

def getCOM(path):
    path_parameter = '/'.join(path.split('/')[:4])+'/parameter/'+path.split('/')[5]+'.dat'
    fp = open(path_parameter)
    dim = fp.readlines()[3].split()
    fp.close()
    temp_positionCOM = []
    with open(path) as f:
        for index, line in enumerate(f):
            if index > 0: temp_positionCOM.append([float(i) for dim_index, i in enumerate(line.split()[1: 4])])
    positionCOM = np.zeros((3, len(temp_positionCOM)))
    for i in range(len(temp_positionCOM)):
        positionCOM[:, i] = temp_positionCOM[i]
    return positionCOM

def getDiffCOM(phi, flow, ncycle, aggregation, nop):
    path1 = '/Users/andrewliu/remote_disk/Data/phi_{}_flow_{}_aggregation_{}KT_ncycle_{}_np_{}/data/sphere_props.0.dat'.format(phi, flow, aggregation, ncycle, nop)
    path2 = '/Users/andrewliu/remote_disk/Data/phi_{}_flow_{}_aggregation_{}KT_ncycle_{}_np_{}/data/sphere_props.1.dat'.format(phi, flow, aggregation, ncycle, nop)
    pos1, pos2 = getCOM(path1), getCOM(path2) 
    diffpos = np.zeros(pos1.shape[1])
    timesteps = list(range(pos1.shape[1]))
    for t in timesteps:
        diffpos[t] = np.linalg.norm(pos1[:, t] - pos2[:, t])
    return(diffpos, timesteps)

for i in range(7):
    try:
        result = getDiffCOM(9.9, i*10+10, 2000, 1, 2)
        plt.subplot(2, 4, i+1)
        plt.plot(result[1], result[0], color = 'b')
        plt.title('flow = {}'.format(i*10+10))
        plt.xlabel('time')
        plt.ylabel('distance')
    except:
        print('no data in flow = {}'.format(i*10+10))
plt.show()
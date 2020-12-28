import numpy as np 
import matplotlib.pyplot as plt
import os
import sys

def getCOM(path):
    temp_positionCOM = []
    with open(path) as f:
        for index, line in enumerate(f):
            if index > 0: temp_positionCOM.append([float(i) for i in line.split()[1: 4]])
    positionCOM = np.zeros((3, len(temp_positionCOM)))
    for i in range(len(temp_positionCOM)):
        positionCOM[:, i] = temp_positionCOM[i]
    return positionCOM

def getDiffCOM(aggregation, flow, ncycle, angle):
    if angle == 0:
        path1 = '/Users/andrewliu/remote_disk/Data/aggregation_{}KT_flow_{}_ncycle_{}_newcutoff/data/sphere_props.0.dat'.format(aggregation, flow, ncycle)
        path2 = '/Users/andrewliu/remote_disk/Data/aggregation_{}KT_flow_{}_ncycle_{}_newcutoff/data/sphere_props.1.dat'.format(aggregation, flow, ncycle)
    else:
        path1 = '/Users/andrewliu/remote_disk/Data/aggregation_{}KT_flow_{}_ncycle_{}_newcutoff_angle_{}/data/sphere_props.0.dat'.format(aggregation, flow, ncycle, angle)
        path2 = '/Users/andrewliu/remote_disk/Data/aggregation_{}KT_flow_{}_ncycle_{}_newcutoff_angle_{}/data/sphere_props.1.dat'.format(aggregation, flow, ncycle, angle)
    pos1, pos2 = getCOM(path1), getCOM(path2) 
    diffpos = np.zeros(pos1.shape[1])
    timesteps = list(range(pos1.shape[1]))
    for t in timesteps:
        diffpos[t] = np.linalg.norm(pos1[:, t] - pos2[:, t])
    return(diffpos, timesteps)

aggregation = sys.argv[1]
flow_list = []
ncycle_list = []
angle_list = []
for fn in os.listdir('/Users/andrewliu/remote_disk/Data'):
    if len(fn.split('_'))>=7:
        if fn.split('_')[1][:-2] == str(aggregation):
            flow_list.append(int(fn.split('_')[3]))
            ncycle_list.append(int(fn.split('_')[5]))
            try: angle_list.append(int(fn.split('_')[-1]))
            except: angle_list.append(0)

sort_index = sorted(range(len(flow_list)), key=lambda k: flow_list[k])
flow_list_sorted = [flow_list[i] for i in sort_index]
ncycle_list_sorted = [ncycle_list[i] for i in sort_index]
angle_list_sorted = [angle_list[i] for i in sort_index]

numcase = len(flow_list)
if numcase > 0:
    column = 6
    if numcase % 6 == 0: row = numcase // column
    else: row = (numcase // column) + 1
else: row, column = 1, numcase

fig = plt.figure(figsize=(36,24))
fig.suptitle('aggregation intensity = {}KT'.format(sys.argv[1]), fontsize = 20)
for i, flow in enumerate(flow_list_sorted):
    ncycle = ncycle_list_sorted[i]
    angle = angle_list_sorted[i]
    try:
        result = getDiffCOM(aggregation, flow, ncycle, angle)
    except:
        print('Warning: No data for aggregation = {}KT, flow = {}, angle = {}'.format(aggregation, flow, angle))   
    plt.subplot(row,column,i+1)
    plt.plot(result[1], result[0], color = 'b')
    plt.title('flow = {}, angle = {}'.format(flow, angle))
    plt.xlabel('time')
    plt.ylabel('distance')
#plt.show()
fig.savefig('/Users/andrewliu/Code/YL_project/Pictures/aggregation_{}KT_newcutoff_angled.png'.format(sys.argv[1]))
print('complete making pic of aggregation = {}KT'.format(sys.argv[1]))
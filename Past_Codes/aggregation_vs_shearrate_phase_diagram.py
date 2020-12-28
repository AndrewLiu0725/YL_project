import numpy as np 
import matplotlib.pyplot as plt
import os
import sys
import time
from scipy import stats

start_time = time.time()

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
    return diffpos


breakups, doublets = [], []
for index, i in enumerate(range(10, 55, 5)):
    print('Fetching data for aggregation = {}KT'.format(i/10.0))
    breakups.append([])
    doublets.append([])

    for fn in os.listdir('/Users/andrewliu/remote_disk/Data'):
        if len(fn.split('_'))>=7:
            if fn.split('_')[1][:-2] == str(i/10.0):
                flow = int(fn.split('_')[3])
                ncycle = int(fn.split('_')[5])
                try: angle = int(fn.split('_')[-1])
                except: angle = 0

                if max(getDiffCOM(i/10.0, flow, ncycle, angle)) > 20: breakups[index].append(flow)
                else: doublets[index].append(flow)
    breakups[index] = list(set(breakups[index]))
    doublets[index] = list(set(doublets[index]))
    doublets[index] = [flow for flow in doublets[index] if flow not in breakups[index]]
    breakups[index].sort()
    doublets[index].sort()
# plot here
fig, ax = plt.subplots()
for marker in ['x', 'o']:
    if marker == 'x':
        x, y = [], []
        for index, i in enumerate(list(range(10, 55, 5))):
            x += list(np.ones(len(breakups[index]))*(i/10.0))
            y += breakups[index]
        ax.scatter(x, y, color = 'r', label = 'break up', marker = marker, s = 100)
    else:
        x, y = [], []
        for index, i in enumerate(list(range(10, 55, 5))):
            x += list(np.ones(len(doublets[index]))*(i/10.0))
            y += doublets[index]
        ax.scatter(x, y, color = 'b', label = 'doublet', marker = marker, s = 100)
m, b, r_value, p_value, std_err = stats.linregress([i/10.0 for i in range(10, 55,5)], [j[0] for j in breakups])
ax.plot(np.linspace(1.0, 5.0, 100), m*np.linspace(1.0, 5.0, 100)+b, label = 'fitted phase transiton', color = 'k', linestyle='dashed', linewidth = 3)
ax.plot([i/10.0 for i in range(10, 55,5)], [j[0] for j in breakups], label = 'phase transition', color = 'g', linewidth = 3)
ax.legend()
plt.xlabel('Aggregation Intensity (KT)', fontsize = 20)
plt.ylabel('Shear Rate'+r'$(\frac {1}{s})$', fontsize = 20)
plt.text(1.0, 140, 'y={}x+{} \n$R^2$={}'.format(m, b, r_value**2), fontsize = 20)
print('Time elpases = ', time.time()-start_time)
plt.show()
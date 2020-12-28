import numpy as np 
import matplotlib.pyplot as plt
import os

stressid = 0
stresss = ['stress_elas_pos', 'stress_inter_v1']
phis = ['504', '610', '706', '811', '987', '1190', '1513']
shearrates = ['937237e-03', '1405855e-02', '1874474e-02', '3748948e-02']

def getstressxy(phi_order, shearrate_order, ratio_order):
    mainpath = '/Users/andrewliu/remote_disk/Data/'
    folderlist = os.listdir(mainpath)
    ratiolist = [] # store the list of ratios for these phi and shear rate
    flag1 = True

    for fn in folderlist:
        fnsplit = fn.split('_')
        try:
            if fnsplit[1] == phis[phi_order]:
                if fnsplit[3] == shearrates[shearrate_order]:
                    ratiolist.append(fnsplit[5])
                    if flag1:
                        temp_fn = fnsplit
                        flag1 = False
        except:
            pass

    # Sort the ratiolist according to their magnitude
    ratiolist = [ratiolist[i] for i in np.argsort(np.array([float(j) for j in ratiolist]))]
    
    temp_fn[5] = ratiolist[ratio_order]
    temp_fn = '_'.join(temp_fn)
    path = mainpath + temp_fn + '/data/{}.dat'.format(stresss[stressid])
    stress_xy = []
    num_timestpes = sum(1 for line in open(path))

    with open(path) as f:  
        for index, line in enumerate(f):
            if index > (num_timestpes - 200):
                stress_xy.append(float(line.split()[1]))
    return stress_xy, ratiolist

stress_xy_rms = np.zeros((3, len(shearrates)))
for k in range(3):
    for j in range(len(shearrates)):
        result = getstressxy(5, j, k)
        stress_xy_rms[k, j] = np.sqrt(np.mean(np.array(result[0])**2))

plt.figure(figsize=(12,10))
plt.xlabel('shear rate', fontsize = 24)
plt.ylabel(r'$\sigma _{xy}$', fontsize = 24)
shearrates_numerical = [float(a[0]+'.'+a[1:]) for a in shearrates]

for i, ratio in enumerate(result[1]):
    plt.plot(shearrates_numerical, stress_xy_rms[i,:], label = ratio[0]+'.'+ratio[1:])
plt.legend(loc = 0)
plt.show()
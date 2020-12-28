import numpy as np 
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import os
import vtk
import time
from vtk import vtkStructuredPointsReader
from vtk.util import numpy_support as VN

starttime = time.time()

def getulsit(volumefraction, shearrate):
        # Get and store the filename of ufile
        # path = '/Users/andrewliu/remote_disk/RBC_omp_single/Data/phi_' + volumefraction + '_rate_' + shearrate + '_np_1/data/'
        path = '/Users/andrewliu/Code/IOP/Data/phi_' + volumefraction + '_rate_' + shearrate + '_np_1/data/'
        datalist = os.listdir(path)
        ulist = []
        timeindex = []
        for fn in datalist:
                if fn[0] == 'u': 
                        ulist.append(fn)
                        timeindex.append(float(fn.split('.')[1]))
        # Sort the ulist accorind to time
        ulist = [ulist[i] for i in np.argsort(np.array(timeindex))]
        return ulist, path

def uofrt(ulist, path, bound, direction):
    ut = np.zeros((len(ulist), 2*bound))
    for index, i in enumerate(ulist):
        reader = vtkStructuredPointsReader()
        reader.SetFileName(path + i)
        reader.ReadAllVectorsOn()
        reader.ReadAllScalarsOn()
        reader.Update()
        data = reader.GetOutput()
        u = VN.vtk_to_numpy(data.GetPointData().GetArray('fluid'))
        dim = list(data.GetDimensions())
        ut[index, :] = u[(dim[0]*dim[1]*int(dim[2]/2)+dim[0]*int(dim[1]/2)): (dim[0]*dim[1]*int(dim[2]/2)+dim[0]*int(dim[1]/2) + 2*bound), direction]
    return ut

bounds = [14, 10]
phis = ['493', '982']
shearrates = ['937237e-03', '1405855e-02', '1874474e-02', '3748948e-02']
shearrate = shearrates[0]
direction = 0 # 0 for x, 1 for y, 2 for z
fieldname = ['ux', 'uy', 'uz']

# Store uz file
for index, phi in enumerate(phis):
    ulist, path = getulsit(phi, shearrate)
    if index == 0:    
        if shearrate == shearrates[0]:
            result1 = uofrt(ulist, path, bounds[index], direction)[1000:,:]
        else:
            result1 = uofrt(ulist, path, bounds[index], direction)[600:,:]
    else:
        result2 = uofrt(ulist, path, bounds[index], direction)[600:,:]
print('Total elapsed time =', time.time() - starttime)

# Make animation here

from matplotlib.animation import FuncAnimation

starttime = time.time()
y_upperbound = max(result1.max(), result2.max())
y_lowerbound = min(result1.min(), result2.min())
fig = plt.figure(figsize = (12, 9))
fig.suptitle('Hydrodynamic Field - {0} \n Shear rate = {1}.{2}'.format(fieldname[direction], shearrate[0], shearrate[1:]), fontsize = 16)
ax = plt.axes(xlim = (-14, 14), ylim = (y_lowerbound, y_upperbound))
line1, line2 = ax.plot([], [],[], [], lw=2)
line1.set_color('darkorange')
line2.set_color('royalblue')

r1 = list(np.arange(-bounds[0], bounds[0]+1))
r1.remove(0)
r2 = list(np.arange(-bounds[1], bounds[1]+1))
r2.remove(0)

def init():
    line1.set_data([],[])
    line2.set_data([],[])
    return tuple([line1, line2])

def animate(i):
    line1.set_data(r1, result1[i, :])
    line2.set_data(r2, result2[i, :])
    ax.set_xlabel('r \n timestep {0}'.format(i), fontsize = 16)
    ax.set_ylabel(fieldname[direction], fontsize = 16)
    ax.legend(['4.93%', '9.82%'])
    return tuple([line1, line2])

anim = FuncAnimation(fig, animate, init_func = init, frames = 200, interval = 100, blit = True)
anim.save('{0}_shearrate_{1}.gif'.format(fieldname[direction], shearrate), writer = 'pillow')
plt.close()
print('Total elapsed time to produce this GIF =', time.time() - starttime)
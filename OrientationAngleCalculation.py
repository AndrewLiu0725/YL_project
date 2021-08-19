# ===============================================================================
# Copyright 2021 An-Jun Liu
# Last Modified Date: 08/04/2021
# ===============================================================================
import numpy as np 
from numpy import linalg as LA
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
from RBC_Utilities import getSuspensionParameterSets
import time
import datetime

"""
This code is to calculate the ensemble averaged orientation angles of suspension system.

The orientation angles are defined as follow:
1. theta is the angle between the major axis of the deformed RBC (minimum principal axis) and the shear direction (x-axis)
2. phi is the angle between the vortex axis (maximum principle axis) and the normal vector (z-axis)

The ensemble average is taken over time steps, particles, and ensemble simulations
"""

'''
To do:
the principal axes is ascending

Candidate for best algo:
data fitting to find best fit plane first and then project all nodes onto the plane and do PCA
'''

bead_number = 642
x_axis = np.array([1, 0, 0])
z_axis = np.array([0, 0, 1])
dim = [144, 24, 144]

# main function
# ===============================================================================
def calcOrientationAngle(path, numberParticle, return_type):
    '''
    return_type: 0 for orientation angles, 1 for principal axes
    '''
    # in the acumulation phase, don't bother to the redundant part
    # just copy the symmetric part after the summation

    # get the data
    data = np.loadtxt(path)
    if return_type:
        principalAxes = np.zeros((numberParticle, 3, 3))
        coms = np.zeros((numberParticle, 3))
    else:
        orientationAngles = np.zeros((numberParticle, 2))

    # run over each partlce
    for i in range(numberParticle):
        inertia_tensor = np.zeros((3, 3))

        # calculate COM
        com = np.sum(data[i*bead_number:(i+1)*bead_number, :], axis = 0)/bead_number

        # run over each node on this particle
        for j in range(bead_number):
            pos_relative = data[i*bead_number+j, :] - com
            inertia_tensor[0, 0] += (pos_relative[1]**2 + pos_relative[2]**2)
            inertia_tensor[0, 1] += (-pos_relative[0]*pos_relative[1])
            inertia_tensor[0, 2] += (-pos_relative[0]*pos_relative[2])
            inertia_tensor[1, 1] += (pos_relative[0]**2 + pos_relative[2]**2)
            inertia_tensor[1, 2] += (-pos_relative[1]*pos_relative[2])
            inertia_tensor[2, 2] += (pos_relative[0]**2 + pos_relative[1]**2)

        # symmetric
        inertia_tensor[1, 0] = inertia_tensor[0, 1]
        inertia_tensor[2, 0] = inertia_tensor[0, 2]
        inertia_tensor[2, 1] = inertia_tensor[1, 2]

        # calculate principle axis
        # the eigenvector is already normalized
        w, v = LA.eig(inertia_tensor)

        if return_type:
            coms[i, :] = [com[j]%dim[j] for j in range(3)]
            principalAxes[i, :, :] = v[np.argsort(v)]
        else:
            orientationAngles[i, 0] = np.arccos(np.dot(v[np.argmax(w)], z_axis)) # phi
            orientationAngles[i, 1] = np.arccos(np.dot(v[np.argmin(w)], x_axis)) # theta

    return [coms, principalAxes] if return_type else orientationAngles



def calcEnsembleAvergaedOrientationAngle(phi, Ca, ensemble_id):
    # get the particle numbers and time steps
    job_name = "h24phi{}Re0.1Ca{}WCA1zero0.8-{}".format(phi, Ca, ensemble_id)
    with open("/userdata4/ajliu/Data_Transfer/" + "{}_parameter.txt".format(job_name)) as f:
        pre_parameters = f.readlines()
    timesteps = int((pre_parameters[pre_parameters.index("timesteps\n")+1])[:-1]) # number of COM files
    particle_numbers = int((pre_parameters[pre_parameters.index("particle_numbers\n")+1])[:-1])

    # get phi and theta time series for each particle
    folder_name = '/raid6/ctliao/Data/HI_ordering/h24_phi{}_Re0.1_Ca{}_WCA1_zero0.8-{}/data/'.format(phi, Ca, ensemble_id)
    st = int(0.9*timesteps)
    orientation_angles = np.zeros((particle_numbers, timesteps-st, 2))
    for t in range(timesteps-st):
        orientation_angles[:, t, :] = calcOrientationAngle(folder_name + 'nodePositions{}.dat'.format(4000*(t+st)), particle_numbers, 0)

    # calculate the ensemble average of the orietation angles (average over time and particles)
    return np.mean(orientation_angles, axis=(0, 1)) # [phi, theta]



def run():
    # function to get ensemble aberaged orientation angle for each set of Ca and volume fraction
    start_time = time.time()

    f_main = open("main.log", "w")

    phi_range = [2.3994, 2.9993, 3.4492, 3.8991, 4.9488, 5.9986]
    #Ca_range = [0.01, 0.02, 0.03, 0.06, 0.07, 0.08, 0.09, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2]
    Ca_range = [0.03, 0.08, 0.1, 0.14, 0.18]

    _, parameter_set = getSuspensionParameterSets()

    total_orientation_angles = np.zeros((len(Ca_range), len(phi_range), 2))
    for Ca_index, Ca in enumerate(Ca_range):
        if Ca == 0.1: time.sleep(10)
        for phi_index, phi in enumerate(phi_range):
            ensemble_count = 0
            for ensemble_id in parameter_set[phi][Ca]:
                try:
                    total_orientation_angles[Ca_index, phi_index, :] += calcEnsembleAvergaedOrientationAngle(phi, Ca, ensemble_id)
                    ensemble_count += 1
                except:
                    pass
            total_orientation_angles[Ca_index, phi_index, :] /= ensemble_count
            f_main.write('(Ca, phi) = ({}, {}) at {}\n'.format(Ca, phi, datetime.datetime.now()))
            f_main.flush()

    np.save('total_orientation_angles.npy', total_orientation_angles)

    f_main.write("Finish!\n")
    f_main.write("Total time elapsed = {}\n".format(str(datetime.timedelta(seconds=time.time()-start_time))))
    f_main.close()



def makePlot(save, filepath, Ca_range, suffix):
    phi_range = np.array([2.3994, 2.9993, 3.4492, 3.8991, 4.9488, 5.9986])
    Ca_range_selected = [0.03, 0.08, 0.1, 0.14, 0.18]

    data = np.load(filepath)

    fig, ax = plt.subplots(figsize = (9, 6))

    for Ca_index, Ca in enumerate(Ca_range):
        if Ca not in Ca_range_selected: continue
        p = ax.plot(phi_range*0.01, data[Ca_index, :, 0]/np.pi, linestyle = '--', marker = "^", label = 'Ca={}, {}={}'.format(Ca, r'$\xi$', r'$\Psi$'))
        ax.plot(phi_range*0.01, data[Ca_index, :, 1]/np.pi, linestyle = '--',marker = "v", color = p[0].get_color(), label = 'Ca = {}, {}={}'.format(Ca, r'$\xi$', r'$\theta$'))
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    ax.set_xlabel(r'$\phi$', fontsize = 12)
    ax.set_ylabel(r'$\dfrac{\langle \xi \rangle }{\pi }$', fontsize = 12)
    fig.tight_layout()
    if save:
        plt.savefig('Pictures/OrientationAngles/OrientationAngle_vs_phi_{}.png'.format(suffix), dpi = 200)
    else:
        plt.show()

    fig, axs = plt.subplots(1, 2, figsize = (15, 6))

    for Ca_index, Ca in enumerate(Ca_range):
        if Ca not in Ca_range_selected: continue
        axs[0].plot(phi_range*0.01, data[Ca_index, :, 0]/np.pi, linestyle = '--', marker = "^", label = 'Ca={}'.format(Ca))
        axs[1].plot(phi_range*0.01, data[Ca_index, :, 1]/np.pi, linestyle = '--',marker = "v", label = 'Ca={}'.format(Ca))
    for i in range(2):
        axs[i].legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        axs[i].set_xlabel(r'$\phi$', fontsize = 12)
        axs[i].set_ylabel(r'$\dfrac{\langle \theta \rangle }{\pi }$' if i else r'$\dfrac{\langle \Psi \rangle }{\pi }$', fontsize = 12)
        axs[i].set_title('{} vs {}'.format(r'$\dfrac{\langle \theta \rangle }{\pi }$' if i else r'$\dfrac{\langle \Psi \rangle }{\pi }$', r'$\phi$'))
    fig.tight_layout()
    if save:
        plt.savefig('Pictures/OrientationAngles/OrientationAngle_vs_phi_seperated_{}.png'.format(suffix), dpi = 200)
    else:
        plt.show()



scale = 10
def checkTwoCell(filepath):
    with open(filepath, 'r') as f:
        data = f.readlines()

    # write corresponding nodePositions file
    numParticle = int(int(data[5].split()[1])/bead_number)
    with open('Data/nodePositions_test.dat', 'w') as f:
        f.writelines(data[6:6+numParticle*bead_number])

    # calculate principal axes and write
    principalAxes = calcOrientationAngle('Data/nodePositions_test.dat', numParticle, 1)
    with open('Data/principalAxes_test.vtk', 'w') as f:
        f.writelines(data[:5])
        f.write('POINTS {} float\n'.format(4*numParticle))
        for i in range(numParticle):
            f.write(" ".join(str(x) for x in (principalAxes[0][i, :] + scale*principalAxes[1][i, 0, :])) + '\n')
            f.write(" ".join(str(x) for x in (principalAxes[0][i, :] - scale*principalAxes[1][i, 0, :])) + '\n')
            f.write(" ".join(str(x) for x in (principalAxes[0][i, :] + scale*principalAxes[1][i, 1, :])) + '\n')
            f.write(" ".join(str(x) for x in (principalAxes[0][i, :] - scale*principalAxes[1][i, 1, :])) + '\n')
        f.write('CELLS {} {}\n'.format(2*numParticle, 3*2*numParticle))
        for i in range(numParticle):
            f.write("2 {} {}\n".format(4*i, 4*i+1))
            f.write("2 {} {}\n".format(4*i+2, 4*i+3))
        f.write('CELL_TYPES {}\n'.format(2*numParticle))
        for i in range(2*numParticle):
            f.write("3\n")

    # remind the user to verify on paraview
    print('Finish!')
    return

def checkSuspension(filepath_vtk, filepath_dat):
    with open(filepath_vtk, 'r') as f:
        data = f.readlines()
    numParticle = int(int(data[5].split()[1])/bead_number)

    # calculate principal axes and write
    principalAxes = calcOrientationAngle(filepath_dat, numParticle, 1)
    with open('Data/principalAxes_test.vtk', 'w') as f:
        f.writelines(data[:5])
        f.write('POINTS {} float\n'.format(4*numParticle))
        for i in range(numParticle):
            f.write(" ".join(str(x) for x in (principalAxes[0][i, :] + scale*principalAxes[1][i, 0, :])) + '\n')
            f.write(" ".join(str(x) for x in (principalAxes[0][i, :] - scale*principalAxes[1][i, 0, :])) + '\n')
            f.write(" ".join(str(x) for x in (principalAxes[0][i, :] + scale*principalAxes[1][i, 1, :])) + '\n')
            f.write(" ".join(str(x) for x in (principalAxes[0][i, :] - scale*principalAxes[1][i, 1, :])) + '\n')
        f.write('CELLS {} {}\n'.format(2*numParticle, 3*2*numParticle))
        for i in range(numParticle):
            f.write("2 {} {}\n".format(4*i, 4*i+1))
            f.write("2 {} {}\n".format(4*i+2, 4*i+3))
        f.write('CELL_TYPES {}\n'.format(2*numParticle))
        for i in range(2*numParticle):
            f.write("3\n")

    # remind the user to verify on paraview
    print('Finish!')
    return

def Axes_PA(data, particleID):
    i = particleID
    inertia_tensor = np.zeros((3, 3))

    # calculate COM
    com = np.sum(data[i*bead_number:(i+1)*bead_number, :], axis = 0)/bead_number

    # run over each node on this particle
    for j in range(bead_number):
        pos_relative = data[i*bead_number+j, :] - com
        inertia_tensor[0, 0] += (pos_relative[1]**2 + pos_relative[2]**2)
        inertia_tensor[0, 1] += (-pos_relative[0]*pos_relative[1])
        inertia_tensor[0, 2] += (-pos_relative[0]*pos_relative[2])
        inertia_tensor[1, 1] += (pos_relative[0]**2 + pos_relative[2]**2)
        inertia_tensor[1, 2] += (-pos_relative[1]*pos_relative[2])
        inertia_tensor[2, 2] += (pos_relative[0]**2 + pos_relative[1]**2)

    # symmetric
    inertia_tensor[1, 0] = inertia_tensor[0, 1]
    inertia_tensor[2, 0] = inertia_tensor[0, 2]
    inertia_tensor[2, 1] = inertia_tensor[1, 2]

    # calculate principle axis
    # the eigenvector is already normalized
    w, v = LA.eig(inertia_tensor)

    return v

def Axes_DF(data, particleID):
    '''
    ax+by+cz=1, (a, b, c) is the normal vector of the fitting plane, 
    which is a candidate for vertex axis
    '''
    i = particleID
    A = np.matrix(data[i*bead_number:(i+1)*bead_number, :])
    b = np.matrix(np.ones((bead_number, 1)))
    fit = (A.T * A).I * A.T * b
    fit = np.array(fit).flatten()
    fit = fit / np.linalg.norm(fit)

    return [fit]

def Axes_DF_Linear(data, particleID):
    '''
    z = ax+by+c, (a, b, 1) is the normal vector of the fitting plane, 
    which is a candidate for vertex axis
    Ax = b --> x = (a, b, c)
    '''
    i = particleID
    x = data[i*bead_number:(i+1)*bead_number, 0]
    y = data[i*bead_number:(i+1)*bead_number, 1]
    z = data[i*bead_number:(i+1)*bead_number, 2]
    A = np.zeros((3, 3))
    A[0, 0] = np.dot(x, x)
    A[0, 1] = np.dot(x, y)
    A[0, 2] = np.sum(x)
    A[1, 0] = A[0, 1]
    A[1, 1] = np.dot(y, y)
    A[1, 2] = np.sum(y)
    A[2, 0] = A[0, 2]
    A[2, 1] = A[1, 2]
    A[2, 2] = bead_number

    b = np.zeros((3, 1))
    b[0, 0] = np.dot(x, z)
    b[1, 0] = np.dot(y, z)
    b[2, 0] = np.sum(z)

    A = np.matrix(A)
    b = np.matrix(b)

    fit = np.matmul(A.I, b)
    fit = np.array(fit).flatten()
    fit[2] = 1
    fit = fit / np.linalg.norm(fit)

    return [fit]


def Axes_PCA(data, particleID):
    i = particleID
    pca = PCA(n_components=3)
    pca.fit(data[i*bead_number:(i+1)*bead_number, :])

    return pca.components_

def test(path, particleID):
    # get the data
    data = np.loadtxt(path)

    # assign paricle index
    i = particleID

    # calculate COM
    com = np.sum(data[i*bead_number:(i+1)*bead_number, :], axis = 0)/bead_number

    # plot the particle
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    ax.scatter(data[i*bead_number:(i+1)*bead_number, 0], data[i*bead_number:(i+1)*bead_number, 1], data[i*bead_number:(i+1)*bead_number, 2])

    c = ['r', 'b', 'g', 'k']
    f = [Axes_PA, Axes_PCA, Axes_DF, Axes_DF_Linear]
    f_n = ['Principal Axes', 'PCA', 'Data Fitting', 'Data Fitting Linear']

    for i in range(4):
        result = f[i](data, i)
        for vector in result:
            ax.plot([com[0]+scale*vector[0], com[0]-scale*vector[0]], 
            [com[1]+scale*vector[1], com[1]-scale*vector[1]], 
            [com[2]+scale*vector[2], com[2]-scale*vector[2]], 
            linewidth = 3, color = c[i])
    #ax.legend()
    plt.show()

### main code
# ===============================================================================
'''
#checkTwoCell('Data/bond0_t1467600.vtk')
#checkSuspension('Data/bond0_t4000000.vtk', 'Data/nodePositions4000000.dat')
Ca_range_1 = [0.01, 0.02, 0.03, 0.06, 0.07, 0.08, 0.09, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2]
Ca_range_2 = [0.03, 0.08, 0.1, 0.14, 0.18]
makePlot(1, 'Data/total_orientation_angles_1.npy', Ca_range_1, 'st_40_et_100')
'''
test('Data/nodePositions4000000.dat', 12)
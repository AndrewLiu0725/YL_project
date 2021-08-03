# ===============================================================================
# Copyright 2021 An-Jun Liu
# Last Modified Date: 08/03/2021
# ===============================================================================
import numpy as np 
from numpy import linalg as LA
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

bead_number = 642
x_axis = np.array([1, 0, 0])
z_axis = np.array([0, 0, 1])

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
        principalAxes = np.zeros((numberParticle, 2, 3))
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
            coms[i, :] = com
            principalAxes[i, 0, :] = v[np.argmax(w)]
            principalAxes[i, 1, :] = v[np.argmin(w)]
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

def makePlot():
    phi_range = np.array([2.3994, 2.9993, 3.4492, 3.8991, 4.9488, 5.9986])
    Ca_range = [0.03, 0.08, 0.1, 0.14, 0.18]

    data = np.load('Data/total_orientation_angles.npy')

    fig, ax = plt.subplots(figsize = (9, 6))

    for Ca_index, Ca in enumerate(Ca_range):
        p = ax.plot(phi_range*0.01, data[Ca_index, :, 0]/np.pi, linestyle = '--', marker = "^", label = 'Ca={}, {}={}'.format(Ca, r'$\xi$', r'$\Psi$'))
        ax.plot(phi_range*0.01, data[Ca_index, :, 1]/np.pi, linestyle = '--',marker = "v", color = p[0].get_color(), label = 'Ca = {}, {}={}'.format(Ca, r'$\xi$', r'$\theta$'))
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    ax.set_xlabel(r'$\phi$', fontsize = 12)
    ax.set_ylabel(r'$\dfrac{\langle \xi \rangle }{\pi }$', fontsize = 12)
    fig.tight_layout()
    plt.savefig('Pictures/OrientationAngles/OrientationAngle_vs_phi.png', dpi = 200)
    #plt.show()

def check(filepath):
    # read bond.vtk

    # write as nodePositionfile

    # calculate principal axes

    # write principal axis vector vtk file

    # remind the user to verify on paraview
    pass


### main code
# ===============================================================================
#calcOrientationAngle('Data/nodePositions_twoCell.dat', 2)
makePlot()
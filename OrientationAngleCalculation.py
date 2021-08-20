# ===============================================================================
# Copyright 2021 An-Jun Liu
# Last Modified Date: 08/20/2021
# ===============================================================================
import numpy as np
import matplotlib.pyplot as plt
from RBC_Utilities import getSuspensionParameterSets
import time
import datetime

"""
This code is to calculate the ensemble averaged orientation angles of suspension system.

The orientation angles are defined as follow:
1. theta is the angle between the major axis of the deformed RBC (minimum principal axis) and the shear direction (x-axis)
2. Psi is the angle between the vortex axis (z-axis) and the normal vector of initial concanve point

The ensemble average is taken over time steps, particles, and ensemble simulations
"""

'''
To do:
the principal axes is ascending

Candidate for best algo:
data fitting to find best fit plane first and then project all nodes onto the plane and do PCA
'''

bead_number = 642
dimple_node = 314
dimple_bonds = [321, 313, 49, 312, 223, 103]
x_axis = np.array([1, 0, 0])
z_axis_positive = np.array([0, 0, 1])
z_axis_negative = np.array([0, 0, -1])
dim = [144, 24, 144]
numOrientationAngles = 1

# main function
# ===============================================================================
def calcOrientationAngle(path, numberParticle):
    # get the data
    data = np.loadtxt(path)
    orientationAngles = np.zeros((numberParticle, numOrientationAngles))

    # run over each partlce
    for i in range(numberParticle):
        normalVectors = np.zeros((6, 3))
        concaveNode = data[dimple_node + i*bead_number, :]
        
        for j in range(6):
            n1 = data[dimple_bonds[j%6]+ i*bead_number, :]-concaveNode
            n2 = data[dimple_bonds[(j+1)%6]+ i*bead_number, :]-concaveNode
            tmpNormalVector = np.cross(n2, n1)
            normalVectors[j, :] = tmpNormalVector/np.linalg.norm(tmpNormalVector)
        
        avgNormalVector = np.sum(normalVectors, axis=0)
        normalizedAvgNormalVector = avgNormalVector/np.linalg.norm(avgNormalVector)
        # Psi
        orientationAngles[i, 0] = min(np.arccos(np.dot(normalizedAvgNormalVector, z_axis_positive)), 
        np.arccos(np.dot(normalizedAvgNormalVector, z_axis_negative)))
    return orientationAngles



def calcEnsembleAvergaedOrientationAngle(phi, Ca, ensemble_id):
    '''
    The ensemble averaged orientation angle is defined as average over time steps and over capsules
    '''
    # get the particle numbers and time steps
    job_name = "h24phi{}Re0.1Ca{}WCA1zero0.8-{}".format(phi, Ca, ensemble_id)
    with open("/userdata4/ajliu/Data_Transfer/" + "{}_parameter.txt".format(job_name)) as f:
        pre_parameters = f.readlines()
    timesteps = int((pre_parameters[pre_parameters.index("timesteps\n")+1])[:-1]) # number of COM files
    particle_numbers = int((pre_parameters[pre_parameters.index("particle_numbers\n")+1])[:-1])

    # get Psi time series for each particle
    folder_name = '/raid6/ctliao/Data/HI_ordering/h24_phi{}_Re0.1_Ca{}_WCA1_zero0.8-{}/data/'.format(phi, Ca, ensemble_id)
    st, et = 40, 100
    orientation_angles = np.zeros((particle_numbers, et-st, numOrientationAngles))
    for t in range(et-st):
        orientation_angles[:, t, :] = calcOrientationAngle(folder_name + 'nodePositions{}.dat'.format(4000*(t+st)), particle_numbers)

    # calculate the ensemble average of the orietation angles (average over time and particles)
    return np.mean(orientation_angles, axis=(0, 1)) # Psi



def run():
    # function to get ensemble aberaged orientation angle for each set of Ca and volume fraction
    start_time = time.time()

    f_main = open("main.log", "w")

    phi_range = [2.3994, 2.9993, 3.4492, 3.8991, 4.9488, 5.9986]
    Ca_range = [0.01, 0.02, 0.03, 0.06, 0.07, 0.08, 0.09, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2]
    #Ca_range = [0.03, 0.08, 0.1, 0.14, 0.18]

    _, parameter_set = getSuspensionParameterSets()

    total_orientation_angles = np.zeros((len(Ca_range), len(phi_range), 1))
    for Ca_index, Ca in enumerate(Ca_range):
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
    Ca_range_selected = [0.03, 0.06, 0.08, 0.1, 0.14, 0.18]

    data = np.load(filepath)

    fig, ax = plt.subplots(figsize = (9, 6))

    for Ca_index, Ca in enumerate(Ca_range):
        if Ca not in Ca_range_selected: continue
        p = ax.plot(phi_range*0.01, data[Ca_index, :, 0]/np.pi, linestyle = '--', marker = "^", label = 'Ca={}, {}={}'.format(Ca, r'$\xi$', r'$\Psi$'))
        #ax.plot(phi_range*0.01, data[Ca_index, :, 1]/np.pi, linestyle = '--',marker = "v", color = p[0].get_color(), label = 'Ca = {}, {}={}'.format(Ca, r'$\xi$', r'$\theta$'))
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    ax.set_xlabel(r'$\phi$', fontsize = 12)
    ax.set_ylabel(r'$\dfrac{\langle \xi \rangle }{\pi }$', fontsize = 12)
    fig.tight_layout()
    if save:
        plt.savefig('Pictures/OrientationAngles/OrientationAngle_vs_phi_{}.png'.format(suffix), dpi = 200)
    else:
        plt.show()
    '''
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
    '''



scale = 10
def checkSuspension(filepath_vtk, filepath_dat):
    with open(filepath_vtk, 'r') as f:
        data = f.readlines()
    numParticle = int(int(data[5].split()[1])/bead_number)
    normalVectors = calcOrientationAngle(filepath_dat, numParticle)

    with open('Data/principalAxes_test.vtk', 'w') as f:
        f.writelines(data[:5])

        f.write('POINTS {} float\n'.format(2*numParticle))
        for i in range(numParticle):
            concaveNode = data[dimple_node + i*bead_number + 6]
            f.write(concaveNode)
            f.write(" ".join([str(float(concaveNode.split()[j]) + scale*normalVectors[i, j]) for j in range(3)]) + '\n')
        
        f.write('CELLS {} {}\n'.format(1*numParticle, 3*numParticle))
        for i in range(numParticle):
            f.write("2 {} {}\n".format(2*i, 2*i+1))

        f.write('CELL_TYPES {}\n'.format(numParticle))
        for i in range(1*numParticle):
            f.write("3\n")

    # remind the user to verify on paraview
    print('Finish!')
    return




### main code
# ===============================================================================
'''
#checkTwoCell('Data/bond0_t1467600.vtk')
#checkSuspension('Data/bond0_t4000000.vtk', 'Data/nodePositions4000000.dat')
Ca_range_1 = [0.01, 0.02, 0.03, 0.06, 0.07, 0.08, 0.09, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2]
Ca_range_2 = [0.03, 0.08, 0.1, 0.14, 0.18]
makePlot(1, 'Data/total_orientation_angles_1.npy', Ca_range_1, 'st_40_et_100')
'''
#print(calcOrientationAngle('Data/nodePositions4000000.dat', 20))
#checkSuspension('Data/bond0_t4000000.vtk', 'Data/nodePositions4000000.dat')
#run()
Ca_range_1 = [0.01, 0.02, 0.03, 0.06, 0.07, 0.08, 0.09, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2]
Ca_range_2 = [0.03, 0.08, 0.1, 0.14, 0.18]
makePlot(0, 'Data/total_orientation_angles.npy', Ca_range_1, 'st_40_et_100')
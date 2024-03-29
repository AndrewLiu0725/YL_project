# ===============================================================================
# Copyright 2021 An-Jun Liu
# Last Modified Date: 08/24/2022
# ===============================================================================
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from RBC_Utilities import getSuspensionParameterSets
import time
import datetime

"""
This code is to calculate the ensemble averaged orientation angles of suspension system.

The orientation angles are defined as follow:
1. theta is the angle between the major axis of the deformed RBC and the shear direction (x-axis)
2. Psi is the angle between the vortex axis (z-axis) and the normal vector of initial concanve point

The ensemble average is taken over time steps, particles, and ensemble simulations
"""

bead_number = 642
dimple_node = [314, 326]
dimple_bonds = [[321, 313, 49, 312, 223, 103], [235, 322, 52, 325, 332, 115]]
x_axis = np.array([1, 0, 0])
z_axis_positive = np.array([0, 0, 1])
z_axis_negative = np.array([0, 0, -1])
dim = [144, 24, 144]
numOrientationAngles = 2

# main function
# ===============================================================================
def findMajorAxisNodes(path, numberParticle):
    # get the data
    data = np.loadtxt(path)
    majorAxisIndices = np.zeros(numberParticle).astype(int)

    # run over each partlce
    for i in range(numberParticle):
        com = np.average(data[i*bead_number:(i+1)*bead_number, :], axis=0)

        # init
        max_distance = np.linalg.norm(data[i*bead_number, :] - com)
        farthest_node = i*bead_number

        # run over each node and compare the to-center-distance
        for j in range(1, bead_number):
            node_index = i*bead_number + j
            tmp_distance = np.linalg.norm(data[node_index, :] - com)
            if tmp_distance > max_distance:
                max_distance = tmp_distance
                farthest_node = node_index

        majorAxisIndices[i] = farthest_node

    return majorAxisIndices



def calcOrientationAngle(path, numberParticle, majorAxisIndices):
    # get the data
    data = np.loadtxt(path)
    orientationAngles = np.zeros((numberParticle, numOrientationAngles))

    # run over each partlce
    for i in range(numberParticle):
        tmpPsi = np.zeros(2)
        for j in range(2):
            normalVectors = np.zeros((6, 3))
            concaveNode = data[dimple_node[j] + i*bead_number, :]

            for k in range(6):
                n1 = data[dimple_bonds[j][k%6]+ i*bead_number, :]-concaveNode
                n2 = data[dimple_bonds[j][(k+1)%6]+ i*bead_number, :]-concaveNode
                tmpNormalVector = np.cross(n2, n1)
                normalVectors[k, :] = tmpNormalVector/np.linalg.norm(tmpNormalVector)
            
            avgNormalVector = np.sum(normalVectors, axis=0)
            normalizedAvgNormalVector = avgNormalVector/np.linalg.norm(avgNormalVector)
            tmpPsi[j] = min(np.arccos(np.dot(normalizedAvgNormalVector, z_axis_positive)), 
            np.arccos(np.dot(normalizedAvgNormalVector, z_axis_negative)))
        # Psi
        orientationAngles[i, 0] = np.average(tmpPsi)

        # theta
        com = np.average(data[i*bead_number:(i+1)*bead_number, :], axis=0)
        majorAxisVector = data[majorAxisIndices[i], :] - com
        x, y = majorAxisVector[0], majorAxisVector[1]
        if x > 0:
            theta = np.arccos(x/np.sqrt(x**2+y**2)) * np.sign(y)
        else:
            theta = (np.pi - np.arccos(x/np.sqrt(x**2+y**2))) * np.sign(-y) 
        orientationAngles[i, 1] = theta

    return orientationAngles # [Psi, theta]



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
    folder_name = '/raid6/ctliao/Data/HI_ordering/h24_phi{}_Re0.1_Ca{}_WCA1_zero0.8-{}/data/'.format(phi, Ca, ensemble_id)

    # find major axis' nodes' indices (use the deformed particle at strain 100)
    majorAxisIndices = findMajorAxisNodes(folder_name + 'nodePositions{}.dat'.format(4000*(100)), particle_numbers)

    # get orientation angles time series for each particle
    st, et = int(0.75*timesteps), timesteps
    orientation_angles = np.zeros((particle_numbers, numOrientationAngles))
    for t in range(st, et):
        orientation_angles += calcOrientationAngle(folder_name + 'nodePositions{}.dat'.format(4000*t), particle_numbers, majorAxisIndices)

    # calculate the ensemble average of the orietation angles (average over particles and time)
    return np.mean(orientation_angles, axis=0)/(et-st) # [Psi, theta]



def run():
    # function to get ensemble aberaged orientation angle for each set of Ca and volume fraction
    start_time = time.time()

    f_main = open("main.log", "w")

    phi_range = [2.3994, 2.9993, 3.4492, 3.8991, 4.9488, 5.9986]
    Ca_range = [0.01, 0.02, 0.03, 0.06, 0.07, 0.08, 0.09, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2]

    _, parameter_set = getSuspensionParameterSets()

    total_orientation_angles = np.zeros((len(Ca_range), len(phi_range), numOrientationAngles))
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

CB_color_cycle = ['blue', 'cyan', 'black',
                  'red', 'purple', 'brown']

def makePlot(save, filepath, Ca_range_data, Ca_range_plot, suffix):
    phi_range = np.array([2.3994, 2.9993, 3.4492, 3.8991, 4.9488, 5.9986])
    data = np.load(filepath)
    
    plt.rcParams['font.family'] = 'DeJavu Serif'
    plt.rcParams['font.serif'] = ['Times New Roman']

    fig, ax = plt.subplots(figsize = (9, 6))

    line_id = 0
    for Ca_index, Ca in enumerate(Ca_range_data):
        if Ca not in Ca_range_plot: continue
        p = ax.plot(phi_range*0.01, data[Ca_index, :, 0]/np.pi, color=CB_color_cycle[line_id], linestyle = '--', marker = "^", label = 'Ca={}, {}={}'.format(Ca, r'$\xi$', r'$\Psi$'))
        ax.plot(phi_range*0.01, data[Ca_index, :, 1]/np.pi, linestyle = '--',marker = "v", color = p[0].get_color(), label = 'Ca = {}, {}={}'.format(Ca, r'$\xi$', r'$\theta$'))
        line_id += 1
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    ax.set_xlabel(r'$\phi$', fontsize = 12)
    ax.set_ylabel(r'$\dfrac{\langle \xi \rangle }{\pi }$', fontsize = 12)
    fig.tight_layout()
    if save:
        plt.savefig('Pictures/OrientationAngles/OrientationAngle_vs_phi_{}.png'.format(suffix), dpi = 200)
    else:
        plt.show()
    
    fig, axs = plt.subplots(1, 2, figsize = (18, 6))

    line_id = 0

    for Ca_index, Ca in enumerate(Ca_range_data):
        if Ca not in Ca_range_plot: continue
        axs[0].plot(phi_range*0.01, data[Ca_index, :, 0]/np.pi, color=CB_color_cycle[line_id], linestyle='--', marker="^", label='Ca={}'.format(Ca))
        axs[1].plot(phi_range*0.01, data[Ca_index, :, 1]/np.pi, color=CB_color_cycle[line_id], linestyle='--', marker="v", label='Ca={}'.format(Ca))
        line_id += 1
    for i in range(2):
        axs[i].legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=15)
        axs[i].set_xlabel(r'$\phi$', fontsize = 20)
        axs[i].set_ylabel(r'$\dfrac{\langle \theta \rangle }{\pi }$' if i else r'$\dfrac{\langle \Psi \rangle }{\pi }$', fontsize = 20)
        #axs[i].set_title('{} vs {}'.format(r'$\dfrac{\langle \theta \rangle }{\pi }$' if i else r'$\dfrac{\langle \Psi \rangle }{\pi }$', r'$\phi$'), fontsize=20)
        axs[i].tick_params(which='both', labelsize=15, width=2, length=8, direction='in')
        axs[i].xaxis.set_major_locator(MultipleLocator(0.01))
        axs[i].xaxis.set_minor_locator(MultipleLocator(0.005))
        #axs[i].yaxis.set_major_locator(MaxNLocator(5)) 
        #axs[i].yaxis.set_minor_locator(MaxNLocator(10))
    fig.tight_layout()
    if save:
        plt.savefig('Pictures/OrientationAngles/OrientationAngle_vs_phi_seperated_{}.png'.format(suffix), dpi = 200)
    else:
        plt.show()
    



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
makePlot(1, 'Data/total_orientation_angles_twoOA.npy', Ca_range_1, Ca_range_2, 'LastQuarter')
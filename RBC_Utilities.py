# ===============================================================================
# Copyright 2021 An-Jun Liu
# Last Modified Date: 06/10/2021
# ===============================================================================
import numpy as np
import matplotlib.pyplot as plt
import sys
import time
from scipy import stats
import ctypes
import os
import logging
import traceback

"""
This utility python file provides functions to calculate doublet fraction, interparticle/elastic stress, intrinsic viscosity, 
and relative viscosity for two-cell and suspension system. 
The Multi-Line Docstring for each function is also provided. 
Please refer to these docstrings before using these utility functions.
"""

# on server set flag = 1, on local set flag = 0
# sevrer: /userdata4/ajliu/ <=>  local: /Users/andrewliu/remote_disk/
flag = 1 

if flag == 1:
    path_preprocess = "/userdata4/ajliu/Data_Transfer/"
    path_CT = "/raid6/ctliao/Data/HI_ordering/"
    path_AJ = "/userdata4/ajliu/RBC_doublet/Data/"

else:
    path_preprocess = "/Users/andrewliu/remote_disk/Data_Transfer/"
    path_CT = "/Users/andrewliu/remote_disk2"
    path_AJ = "/Users/andrewliu/remote_disk/RBC_doublet/Data/"


def calcDoubletFraction(input_phi, input_Ca, input_criteria_T, input_criteria_Dms, make_plot, depend, system):
    """
    Input:
    
    input_criteria_Dms must be either lists or floats
    
    make_plot = 0 means doesn't need to make plot, 1 means saving df vs t plot without playing it, 2 means displaying the plot

    system = 0 means two-cell system, 1 means suspension

    depend means angle when system = 0, esemble id when system = 1
    
    Output:
    
    [[doublet fraction], timesteps, [state series], diffpos, particle_numbers]
    
    In state series, 1 means doublet, 2 means kayaking, 3 means sliding, and 4 means transition states

    Exit code:
    0 Success
    1 OSError

    Note:
    1. All the arrays shared by python and C program are stored in double so that we don't need to consider the data compatibility.
    2. You have to preprocess the data before using this utility function.
    """

    # C extension
    lib = ctypes.cdll.LoadLibrary('./cforDoublet_Functions.so')

    # Parameter
    # ===============================================================================
    phi = input_phi
    Ca = input_Ca
    Dm = 15.64
    #criteria_Ts = input_criteria_Ts if type(input_criteria_Ts) is list else [input_criteria_Ts]
    criteria_Dms = input_criteria_Dms if type(input_criteria_Dms) is list else [input_criteria_Dms]
    
    # Two-cell system
    if system == 0:
        angle = depend
        ncycle = 4000
        job_name = "phi_{}_Re_0.1_Ca_{}_aggregation_1KT_ncycle_{}_np_2_angle_{}".format(phi, Ca, ncycle, angle)
    
    # Suspension system
    else:
        ensemble_id = depend
        job_name = "h24phi{}Re0.1Ca{}WCA1zero0.8-{}".format(phi, Ca, ensemble_id)

    
    # Read the preprocessed files
    # ===============================================================================
    try:
        # Read the parameters
        with open(path_preprocess + "{}_parameter.txt".format(job_name)) as f:
            pre_parameters = f.readlines()
        timesteps = int((pre_parameters[pre_parameters.index("timesteps\n")+1])[:-1]) # number of COM files
        particle_numbers = int((pre_parameters[pre_parameters.index("particle_numbers\n")+1])[:-1])
        points_per_particle = int((pre_parameters[pre_parameters.index("points_per_particle\n")+1])[:-1])
        interval = int((pre_parameters[pre_parameters.index("interval\n")+1])[:-1]) # number of Ypos_t files
        
        dim = np.zeros(3, dtype = np.int32)
        tmp_dim = eval((pre_parameters[pre_parameters.index("dim\n")+1])[:-1])
        for i in range(3):
            dim[i] = tmp_dim[i]

        # Read the position of center of mass
        COMs = np.load(path_preprocess + "{}_COMs.npy".format(job_name))
        COMs_NB = np.load(path_preprocess + "{}_COMs_NB.npy".format(job_name))

        # Read the y positions
        Ypos_t = np.load(path_preprocess + "{}_Ypos_t.npy".format(job_name))
    
    except OSError:
        logging.error(" calcDoubletFraction():\nNo preprocessed data for simulation: {}\n".format(job_name))
        sys.exit(1)


    # Note: All the following computations can be parallelized

    # Calculate t_rot
    # ===============================================================================
    # Format of Ypos_t is Ypos_t[t, node_id]
    Periods = np.zeros(particle_numbers*points_per_particle)
    offset = int(0.0125*timesteps) # igonre frequency lower than 0.0125, i.e. rotation time > 80 starins

    for i in range(particle_numbers*points_per_particle):
        P = np.fft.rfft(Ypos_t[:, i])
        Periods[i] = (interval)/(np.argmax(np.abs(P[offset:]))+offset) # divide timesteps by 2 since here is the number of the timesteps of bond0.vtk

    # convert from the time unit in nodeposition format to COM format
    # in old output format, time unit in nodeposition format = WriteConfig and time unit in COM format = WriteProps
    # in new output format, time unit in nodeposition format = time unit in COM format = WriteProps
    # Updated Date: 01/29/2021 by An-Jun Liu
    #rotation_time = stats.gmean(Periods)*(timesteps/interval) 
    rotation_time = stats.trim_mean(Periods, 0.1)*(timesteps/interval) # remove the possible outliers



    # Calculate pair distance
    # ===============================================================================
    # Format of COMs is COMs[particle_id, t, 3]
    
    number_of_pairs = int((particle_numbers-1)*particle_numbers/2)
    diffpos = np.zeros((number_of_pairs, timesteps), dtype = np.float64)
    uncorrected_diffpos = np.zeros((number_of_pairs, timesteps), dtype = np.float64)
    indice_pairs = np.zeros((number_of_pairs, 2), dtype = np.int32)
    
    count = 0
    for i in range(particle_numbers-1):
        for j in range(i+1, particle_numbers):
            indice_pairs[count, 0], indice_pairs[count, 1] = i, j
            uncorrected_diffpos[count, :] = np.linalg.norm((COMs[i, :, :] - COMs[j, :, :]), axis=1)
            count += 1

    c_correctDiffpos = lib.correctDiffpos
    c_correctDiffpos(ctypes.c_void_p(diffpos.ctypes.data), ctypes.c_void_p(COMs.ctypes.data), ctypes.c_int(number_of_pairs),
    ctypes.c_int(timesteps), ctypes.c_double(Dm), ctypes.c_void_p(indice_pairs.ctypes.data), ctypes.c_void_p(dim.ctypes.data)) 
    

          
    # Calculate doublet fraction
    # ===============================================================================
    if (make_plot == 1  or make_plot == 2):
        fig, ax1 = plt.subplots(figsize = (8, 6))

    output_DF, output_state = [], []
    criteria_T = input_criteria_T
    period = int(round(criteria_T*rotation_time))
    #print("period =", period)
    for criteria_Dm in criteria_Dms:
        doublet_or_not = np.zeros((number_of_pairs, timesteps - period), dtype = np.int32) # 1 means there is a doublet
        state_series = np.zeros((number_of_pairs, timesteps - period), dtype = np.int32)
        end_time = np.zeros(1, dtype = np.int32)

        c_calcDF = lib.calcDF

        c_calcDF(ctypes.c_void_p(doublet_or_not.ctypes.data), ctypes.c_void_p(diffpos.ctypes.data), ctypes.c_int(period),
        ctypes.c_int(timesteps), ctypes.c_int(number_of_pairs), ctypes.c_double(Dm), ctypes.c_double(criteria_Dm),
        ctypes.c_void_p(COMs_NB.ctypes.data), ctypes.c_void_p(indice_pairs.ctypes.data), ctypes.c_void_p(state_series.ctypes.data),
        ctypes.c_void_p(end_time.ctypes.data), ctypes.c_void_p(dim.ctypes.data))
        
        # Assume that there is no aggregation structure which consists more than 2 partilces
        number_of_doublets = np.sum(doublet_or_not, axis = 0)
        # one doublet has two RBCs so to calculate doublet fraction, we need multiply # of doublets by 2
        output_DF.append(number_of_doublets*2/particle_numbers)
        output_state.append(state_series)
        
        if (make_plot == 1  or make_plot == 2):
            ax1.plot(np.array(list(range(timesteps - period)))/2, number_of_doublets*2/particle_numbers, label = "{}Dm, {}t_rot".format(criteria_Dm, criteria_T))
        


    # Make plot
    # ===============================================================================
    if (make_plot == 1 or make_plot == 2):
        ax1.set_xlabel("timesteps (in bond.vtk unit)", fontsize = 20)
        ax1.set_ylabel("doublet fraction", fontsize = 20)
        if system == 0:
            ax1.set_title("Two-cell system\nphi={}, Ca={}, angle={}, Re=0.1, h=24".format(phi, Ca, angle), fontsize = 20)
        else:
            ax1.set_title("Suspension\nphi={}, Ca={}, esemble {}\nRe=0.1, D=1, eqWCA=0.8, h=24".format(phi, Ca, ensemble_id), fontsize = 20)
        ax1.legend(fontsize = 16, bbox_to_anchor=(1.1,1.2))
        ax2 = ax1.twinx()
        mn, mx = ax1.get_ylim()
        ax2.set_ylim(mn*particle_numbers/2, mx*particle_numbers/2)
        ax2.set_ylabel('# of doublets', fontsize = 20)
        fig.tight_layout()
        if make_plot == 1:
            plt.savefig("./Pictures/{}_DoubletFractionTimeSeries.png".format(job_name), dpi = 300)
            plt.close()
        else:
            plt.show()

    #return [output_DF, timesteps, output_state, diffpos, COMs, COMs_NB]
    return [output_DF, end_time[0], output_state, diffpos, particle_numbers]



def getStress(input_phi, input_Ca, stress_category_id, ncycle, depend, system):
    """
    Input:

    input_phi is in percentage unit

    stress_category_id = 0 means elastic stress tensor, 1 means inter-particle stress tensor

    ncycle is only meaningful for two-cell system, i.e. system = 1

    system = 0 means two-cell system, 1 means suspension

    depend means angle when system = 0, esemble id when system = 1

    Output:

    An ndarray of shape (t, 9), i.e. time series for each component

    Exit code:
    0 Success
    1 StopIteration
    2 OSError
    3 Others
    """

    # Parameters
    # ===============================================================================
    stress_category = ['stress_elas_pos', 'stress_inter']
    phi = input_phi
    Ca = input_Ca

    # Get file path
    # ===============================================================================
    # Two-cell system
    if system == 0:
        angle = depend
        path = path_AJ + "phi_{}_Re_0.1_Ca_{}_aggregation_1KT_ncycle_{}_np_2_angle_{}/data/{}.dat".format(phi, Ca, ncycle, angle, stress_category[stress_category_id])
    
    # Suspesion system
    else:
        ensemble_id = depend
        path = path_CT + "h24_phi{}_Re0.1_Ca{}_WCA1_zero0.8-{}/data/{}.dat".format(phi, Ca, ensemble_id, stress_category[stress_category_id])
        if not os.path.isfile(path):
            path = path_CT + "h24phi{}Re0.1Ca{}WCA1zero0.8-{}/data/{}.dat".format(phi, Ca, ensemble_id, stress_category[stress_category_id])
    

    # Output and Error Handling
    # ===============================================================================
    try:
        stress = np.loadtxt(path, skiprows = 2)

    except StopIteration:
        logging.error(" getStress():\nWrong value of timestep_start/timestep_end: {}\n".format(path))
        sys.exit(1)

    except OSError:
        logging.error(" getStress():\nNo such file: {}\n".format(path))
        sys.exit(2)

    except Exception:
        logging.error(traceback.format_exc())
        sys.exit(3)

    return stress



def getIntrinsicViscosity(input_phi, input_Ca, ncycle, depend, system):
    """
    Input:

    input_phi is in percentage unit

    ncycle is only meaningful for two-cell system, i.e. system = 1

    system = 0 means two-cell system, 1 means suspension

    depend means angle when system = 0, esemble id when system = 1

    Output:

    An numpy array of length t, i.e. the time series of [eta]

    Exit code:
    0 Success
    1 StopIteration
    2 OSError
    3 Others
    """

    # Parameters
    # ===============================================================================
    eta_f = 6.0
    phi = input_phi
    Ca = input_Ca

    # Get file path
    # ===============================================================================
    # Two-cell system
    if system == 0:
        angle = depend
        path = path_AJ + "phi_{}_Re_0.1_Ca_{}_aggregation_1KT_ncycle_{}_np_2_angle_{}/data/wallStress.dat".format(phi, Ca, ncycle, angle)
    # Suspension system
    else:
        ensemble_id = depend
        path = path_CT + "h24_phi{}_Re0.1_Ca{}_WCA1_zero0.8-{}/data/wallStress.dat".format(phi, Ca, ensemble_id)
        if not os.path.isfile(path):
            path = path_CT + "h24phi{}Re0.1Ca{}WCA1zero0.8-{}/data/wallStress.dat".format(phi, Ca, ensemble_id)

    # Output and Error Handling
    # ===============================================================================
    try:
        #eta = np.loadtxt(path, skiprows = int(timestep_start+1), max_rows = int(timestep_end-timestep_start))[:,1]
        eta = np.loadtxt(path, skiprows = 1)[:,1]
        intrinsic_eta = (eta-eta_f)/(eta_f*input_phi*0.01) # phi is in percentage unit

    except StopIteration:
        logging.error(" getIntrinsicViscosity():\nWrong value of timestep_start/timestep_end: {}\n".format(path))
        sys.exit(1)

    except OSError:
        logging.error(" getIntrinsicViscosity():\nNo such file: {}\n".format(path))
        sys.exit(2)

    except Exception:
        logging.error(traceback.format_exc())
        sys.exit(3)

    return intrinsic_eta



def getRelativeViscosity(input_phi, input_Ca, ncycle, depend, system):
    """
    Input:

    input_phi is in percentage unit

    ncycle is only meaningful for two-cell system, i.e. system = 1

    system = 0 means two-cell system, 1 means suspension

    depend means angle when system = 0, esemble id when system = 1

    Output:
    An numpy array of length t, i.e. the time series of eta/eta_f

    Exit code:
    0 Success
    1 StopIteration
    2 OSError
    3 Others
    """

    # Parameters
    # ===============================================================================
    eta_f = 6.0
    phi = input_phi
    Ca = input_Ca

    # Get file path
    # ===============================================================================
    # Two-cell system
    if system == 0:
        angle = depend
        path = path_AJ + "phi_{}_Re_0.1_Ca_{}_aggregation_1KT_ncycle_{}_np_2_angle_{}/data/wallStress.dat".format(phi, Ca, ncycle, angle)
    # Suspension system
    else:
        ensemble_id = depend
        path = path_CT + "h24_phi{}_Re0.1_Ca{}_WCA1_zero0.8-{}/data/wallStress.dat".format(phi, Ca, ensemble_id)
        if not os.path.isfile(path):
            path = path_CT + "h24phi{}Re0.1Ca{}WCA1zero0.8-{}/data/wallStress.dat".format(phi, Ca, ensemble_id)
    
    # Output and Error Handling
    # ===============================================================================
    try:
        #eta = np.loadtxt(path, skiprows = int(timestep_start+1), max_rows = int(timestep_end-timestep_start))[:,1]
        eta = np.loadtxt(path, skiprows = 1)[:,1]
        relative_eta = eta/eta_f

    except StopIteration:
        logging.error(" getRelativeViscosity():\nWrong value of timestep_start/timestep_end: {}\n".format(path))
        sys.exit(1)

    except OSError:
        logging.error(" getRelativeViscosity():\nNo such file: {}\n".format(path))
        sys.exit(2)

    except Exception:
        logging.error(traceback.format_exc())
        sys.exit(3)

    return relative_eta



def parser(string):
    """
    This parser is utility function for getSuspensionParameterSets().
    """
    # output [phi, Ca, ensemble id]
    string = string.split("-")
    ensemble_id = int(string[1])
    # format: h24_phi4.4989_Re0.1_Ca0.06_WCA1_zero0.8-8
    if "_" in string[0]:
        string = string[0].split("_")
        phi = float(string[1][3:])
        Ca = float(string[3][2:])
    # format: h24phi4.9488Re0.1Ca0.06WCA1zero0.8-4
    else:
        phi = float(string[0][string[0].find('i')+1:string[0].find('R')])
        Ca = float(string[0][string[0].find('a')+1:string[0].find('W')])
    return [phi, Ca, ensemble_id]



def getSuspensionParameterSets():
    """
    Output:

    [phis, parameter_set]

    phis: a list of phis which will be used

    parameter_set: a two-layered dictionary which stores the sets of parameters.
    
    Usage: parameter_set[phi][Ca] = [ensemble_id]
    """
    path = "/raid6/ctliao/Data/HI_ordering/"

    parameter_set = {}
    # create parameter_set (two layer dict)
    for fn in os.listdir(path):
        if (fn[0] == "h") and (os.path.isdir(path+fn)):
            result = parser(fn)
            [phi, Ca, ensemble_id] = result
            if phi in parameter_set.keys():
                if Ca in parameter_set[phi].keys():
                    parameter_set[phi][Ca].append(ensemble_id)
                else:
                    parameter_set[phi][Ca] = [ensemble_id]
            else:
                parameter_set[phi] = {}
                parameter_set[phi][Ca] = [ensemble_id]


    phis = []

    for phi in parameter_set.keys():
        if (len(parameter_set[phi].keys()) > 5):
            phis.append(phi)

    return [phis, parameter_set]
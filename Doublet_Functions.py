import numpy as np 
import matplotlib.pyplot as plt
import sys
import time
from scipy import stats
import ctypes

def calcDoubletFraction(input_phi, input_Ca, input_criteria_T, input_criteria_Dms, make_plot, *input_others):
    """
    Input:
    
    input_criteria_Dms must be either lists or floats
    
    make_plot = 1 means making df vs t plot without displaying it, 2 means displaying it, 0 means doesn't need to make plot
    
    Output:
    
    [[doublet fraction], timesteps, [state series], diffpos, particle_numbers]
    
    In state series, 1 means doublet, 2 means kayaking, 3 means sliding, and 4 means transition states

    Note: All the arrays shared by python and C program are stored in double so that we don't need to consider the data compatibility
    """

    # C extension
    lib = ctypes.cdll.LoadLibrary('./cforDoublet_Functions.so')

    # Parameter
    ###########################################################################
    phi = input_phi
    Ca = input_Ca
    Dm = 15.64
    #criteria_Ts = input_criteria_Ts if type(input_criteria_Ts) is list else [input_criteria_Ts]
    criteria_Dms = input_criteria_Dms if type(input_criteria_Dms) is list else [input_criteria_Dms]
    
    # Two-cell system
    if len(input_others) > 0:
        WriteProps = 3669
        angle = input_others[0]
        ncycle = 2000
        job_name = "phi_{}_Re_0.1_Ca_{}_aggregation_1KT_ncycle_{}_np_2_angle_{}".format(phi, Ca, ncycle, angle)
    
    # Suspension system
    else:
        job_name = "h24phi{}Re0.1Ca{}WCA1zero0.8".format(phi, Ca)

    
    # Read the preprocessed files by mounting the folder in server, use eval() function
    ###########################################################################
    # Read the parameters
    with open("/Users/andrewliu/remote_disk/Data_Transfer/{}_parameter.txt".format(job_name)) as f:
        pre_parameters = f.readlines()
    timesteps = int((pre_parameters[pre_parameters.index("timesteps\n")+1])[:-1]) # number of COM files
    particle_numbers = int((pre_parameters[pre_parameters.index("particle_numbers\n")+1])[:-1])
    points_per_particle = int((pre_parameters[pre_parameters.index("points_per_particle\n")+1])[:-1])
    interval = int((pre_parameters[pre_parameters.index("interval\n")+1])[:-1]) # number of Ypos_t files
    
    dim = np.zeros(3, dtype = np.int32)
    tmp_dim = eval((pre_parameters[pre_parameters.index("dim\n")+1])[:-1])
    for i in range(3):
        dim[i] = tmp_dim[i]

    # Read WriteProps for suspension system
    if len(input_others) == 0:
        WriteProps = int((pre_parameters[pre_parameters.index("WriteProps\n")+1])[:-1])
    
    

    # Read the position of center of mass
    COMs = np.load("/Users/andrewliu/remote_disk/Data_Transfer/{}_COMs.npy".format(job_name))
    COMs_NB = np.load("/Users/andrewliu/remote_disk/Data_Transfer/{}_COMs_NB.npy".format(job_name))

    # Read the y positions
    Ypos_t = np.load("/Users/andrewliu/remote_disk/Data_Transfer/{}_Ypos_t.npy".format(job_name))


    # All the following computations can be parallelized

    # Calculate t_rot
    ###########################################################################
    # Format of Ypos_t is Ypos_t[t, node_id]
    Periods = np.zeros(particle_numbers*points_per_particle)

    for i in range(particle_numbers*points_per_particle):
        offset = 10 # igonre the physically impossible high frequency regime
        P = np.fft.rfft(Ypos_t[:, i])
        Periods[i] = (interval)/(np.argmax(np.abs(P[offset:]))+offset) # divide timesteps by 2 since here is the number of the timesteps of bond0.vtk

    rotation_time = stats.gmean(Periods)*(timesteps/interval) # in unit of strain, 2 comes from the ratio WriteConfig/WriteProps (= 4000/2000 = 2)



    # Calculate pair distance
    ###########################################################################
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
    ###########################################################################
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
    ###########################################################################
    if (make_plot == 1  or make_plot == 2):
        ax1.set_xlabel("timesteps (in bond.vtk unit)", fontsize = 20)
        ax1.set_ylabel("doublet fraction", fontsize = 20)
        if len(input_others) > 0:
            ax1.set_title("Two-cell system\nphi={}, Ca={}, angle={}, Re=0.1, h=24".format(phi, Ca, angle), fontsize = 20)
        else:
            ax1.set_title("Suspension\nphi={}, Ca={}, Re=0.1, D=1, eqWCA=0.8, h=24".format(phi, Ca), fontsize = 20)
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

def getStress(input_phi, input_Ca, stress_category_id, timestep_start, timestep_end, *input_others):
    """
    Input:
    stress_category_id = 0 means elastic stress tensor, 1 means inter-particle stress tensor
    Output:
    An ndarray of shape (t, 9), i.e. time series for each component
    """
    stress_category = ['stress_elas_pos', 'stress_inter']
    phi = input_phi
    Ca = input_Ca

    # Two-cell system
    if len(input_others) > 0:
        angle = input_others[0]
        path = ("/Users/andrewliu/remote_disk/RBC_doublet/Data/phi_{}_Re_0.1_Ca_{}_aggregation_1KT_ncycle_2000_np_2_angle_{}/data/{}.dat"
                .format(phi, Ca, angle, stress_category[stress_category_id]))
    
    # Suspesion system
    else:
        if phi == 4:
            path = "/Users/andrewliu/remote_disk/phi4/eqWCA0.8/h24phi4Re0.1Ca{}D0eqWCA0.8/data/{}.dat".format(Ca, stress_category[stress_category_id])
        elif phi == 5:
            path = "/Users/andrewliu/remote_disk/phi5/Re0.1/eqWCA0.8/h24phi5Re0.1Ca{}D0eqWCA0.8/data/{}.dat".format(Ca, stress_category[stress_category_id])
    
    stress = np.loadtxt(path, skiprows = timestep_start+2, max_rows = (timestep_end-timestep_start))  
    return stress

def getInstrinsicViscosity(input_phi, input_Ca, timestep_start, timestep_end, *input_others):
    """
    Output:
    An numpy array of length t, i.e. the time series of [eta]
    """
    # Parameters
    eta_f = 6.0
    phi = input_phi
    Ca = input_Ca

    # Two-cell system
    if len(input_others) > 0:
        angle = input_others[0]
        path = ("/Users/andrewliu/remote_disk/RBC_doublet/Data/phi_{}_Re_0.1_Ca_{}_aggregation_1KT_ncycle_2000_np_2_angle_{}/data/wallStress.dat"
                .format(phi, Ca, angle))
    # Suspension system
    else:
        if phi == 4:
            path = "/Users/andrewliu/remote_disk/phi4/eqWCA0.8/h24phi4Re0.1Ca{}D0eqWCA0.8/data/wallStress.dat".format(Ca)
        elif phi == 5:
            path = "/Users/andrewliu/remote_disk/phi5/Re0.1/eqWCA0.8/h24phi5Re0.1Ca{}D0eqWCA0.8/data/wallStress.dat".format(Ca)

    eta = np.loadtxt(path, skiprows = timestep_start+1, max_rows = (timestep_end-timestep_start))[:,1]
    intrinsic_eta = (eta-eta_f)/(eta_f*input_phi*0.01) # phi is in percentage unit
    return intrinsic_eta

def getRelativeViscosity(input_phi, input_Ca, timestep_start, timestep_end, *input_others):
    """
    Output:
    An numpy array of length t, i.e. the time series of eta/eta_f
    """
    # Parameters
    eta_f = 6.0
    phi = input_phi
    Ca = input_Ca

    # Two-cell system
    if len(input_others) > 0:
        angle = input_others[0]
        path = ("/Users/andrewliu/remote_disk/RBC_doublet/Data/phi_{}_Re_0.1_Ca_{}_aggregation_1KT_ncycle_2000_np_2_angle_{}/data/wallStress.dat"
                .format(phi, Ca, angle))
    # Suspension system
    else:
        if phi == 4:
            path = "/Users/andrewliu/remote_disk/phi4/eqWCA0.8/h24phi4Re0.1Ca{}D0eqWCA0.8/data/wallStress.dat".format(Ca)
        elif phi == 5:
            path = "/Users/andrewliu/remote_disk/phi5/Re0.1/eqWCA0.8/h24phi5Re0.1Ca{}D0eqWCA0.8/data/wallStress.dat".format(Ca)

    eta = np.loadtxt(path, skiprows = timestep_start+1, max_rows = (timestep_end-timestep_start))[:,1]
    relative_eta = eta/eta_f
    return relative_eta
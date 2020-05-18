import numpy as np 
import matplotlib.pyplot as plt
import sys
import time
from scipy import stats
import ctypes

def calcDoubletFraction(input_phi, input_Ca, input_criteria_Ts, input_criteria_Dms, make_plot, *input_others):
    """
    Input:
    input_criteria_Ts and input_criteria_Dms must be either list or a float
    make_plot = 1 means making plots, 0 means doesn't need to make plots
    
    Output:
    [[strict doublet fraction], [general doublet fraction], timesteps (COMs)]

    All the arrays shared by python and C program are stored in double so that we don't need to consider the data compatibility
    """

    # C extension
    lib = ctypes.cdll.LoadLibrary('./cforDoublet_Functions.so')

    # Parameter
    ###########################################################################
    phi = input_phi
    Ca = input_Ca
    D = 0
    eqWCA = 0.8
    Dm = 15.64
    criteria_Ts = input_criteria_Ts if type(input_criteria_Ts) is list else [input_criteria_Ts]
    criteria_Dms = input_criteria_Dms if type(input_criteria_Dms) is list else [input_criteria_Dms]
    
    # Two-cell system
    if len(input_others) > 0:
        WriteProps = 3669
        angle = input_others[0]
        ncycle = 2000
        job_name = "phi_{}_Re_0.1_Ca_{}_aggregation_1KT_ncycle_{}_np_2_angle_{}".format(phi, Ca, ncycle, angle)
    
    # Suspension system
    else:
        WriteProps = 2000
        job_name = "h24phi{}Re0.1Ca{}D{}eqWCA{}".format(phi, Ca, D, eqWCA)

    
    # Read the preprocessed files by mounting the folder in server, use eval() function
    ###########################################################################
    # Read the parameters
    with open("/Users/andrewliu/remote_disk/Data_Transfer/{}_parameter.txt".format(job_name)) as f:
        pre_parameters = f.readlines()
    timesteps = int((pre_parameters[pre_parameters.index("timesteps\n")+1])[:-1])
    particle_numbers = int((pre_parameters[pre_parameters.index("particle_numbers\n")+1])[:-1])
    points_per_particle = int((pre_parameters[pre_parameters.index("points_per_particle\n")+1])[:-1])
    
    dim = np.zeros(3, dtype = np.int32)
    # Two-cell system
    try:
        tmp_dim = eval((pre_parameters[pre_parameters.index("dim\n")+1])[:-1])
        for i in range(3):
            dim[i] = tmp_dim[i]
    # Suspension system
    except:
        for i in range(3):
            dim[i] = [144, 24, 144][i]
    

    # Read the position of center of mass
    COMs = np.load("/Users/andrewliu/remote_disk/Data_Transfer/{}_COMs.npy".format(job_name))

    # Read the y positions
    Ypos_t = np.load("/Users/andrewliu/remote_disk/Data_Transfer/{}_Ypos_t.npy".format(job_name))


    # All the following computations can be parallelized

    # Calculate t_rot
    ###########################################################################
    # Format of Ypos_t is Ypos_t[t, node_id]
    Periods = np.zeros(particle_numbers*points_per_particle)

    for i in range(particle_numbers*points_per_particle):
        #Ypos_t_norm = Ypos_t[:, i] - np.mean(Ypos_t[:, i])
        offset = 10
        P = np.fft.rfft(Ypos_t[:, i])
        Periods[i] = (timesteps/2)/(np.argmax(np.abs(P[offset:]))+offset) # divide timesteps by 2 since here is the number of the timesteps of bond0.vtk
  
    rotation_time = stats.gmean(Periods)*2 # in unit of strain, 2 comes from the ratio WriteConfig/WriteProps (= 4000/2000 = 2)
    
    if (make_plot == 1):
        plt.figure(2)
        plt.hist(Periods, bins = particle_numbers)
        plt.xlabel("period")
        plt.ylabel("number")
        plt.title("Geometric mean = {}\nstd = {}".format(rotation_time/2, np.std(Periods)))
        #plt.savefig("./Pictures/Period_histogram_{}.png".format(job_name))
        #print("Rotation time = ",rotation_time)
        plt.close()



    # Calculate pair distance
    ###########################################################################
    # Format of COMs is COMs[particle_id, t, 3]
    
    number_of_pairs = int((particle_numbers-1)*particle_numbers/2)
    diffpos = np.zeros((number_of_pairs, timesteps), dtype = np.float64)
    indice_pairs = np.zeros((number_of_pairs, 2), dtype = np.int32)

    count = 0
    for i in range(particle_numbers-1):
        for j in range(i+1, particle_numbers):
            indice_pairs[count, 0], indice_pairs[count, 1] = i, j
            diffpos[count, :] = np.linalg.norm((COMs[i, :, :] - COMs[j, :, :]), axis=1)
            count += 1       

    # Correct diffpos
    c_correctDiffpos = lib.correctDiffpos
    c_correctDiffpos(ctypes.c_void_p(diffpos.ctypes.data), ctypes.c_void_p(COMs.ctypes.data), ctypes.c_int(number_of_pairs),
    ctypes.c_int(timesteps), ctypes.c_double(Dm), ctypes.c_void_p(indice_pairs.ctypes.data), ctypes.c_void_p(dim.ctypes.data)) 

    
          
    # Calculate doublet fraction
    ###########################################################################
    if (make_plot == 1):
        fig, ax1 = plt.subplots(figsize = (16, 12))

    output_g = []
    output_s = []
    for criteria_Dm in criteria_Dms:
        for criteria_T in criteria_Ts:
            period = int(round(criteria_T*rotation_time))

            doublet_or_not_strict = np.zeros((number_of_pairs, timesteps - period), dtype = np.int32) # 1 means there is a doublet
            doublet_or_not_general = np.zeros((number_of_pairs, timesteps - period), dtype = np.int32) # 1 means there is a doublet
  
            c_calcDFGeneral = lib.calcDFGeneral
            c_calcDFStrict = lib.calcDFStrict

            c_calcDFGeneral(ctypes.c_void_p(doublet_or_not_general.ctypes.data), ctypes.c_void_p(diffpos.ctypes.data), ctypes.c_int(period),
            ctypes.c_int(timesteps), ctypes.c_int(number_of_pairs), ctypes.c_double(Dm), ctypes.c_double(criteria_Dm))

            c_calcDFStrict(ctypes.c_void_p(doublet_or_not_strict.ctypes.data), ctypes.c_void_p(diffpos.ctypes.data), ctypes.c_int(period),
            ctypes.c_int(timesteps), ctypes.c_int(number_of_pairs), ctypes.c_double(Dm), ctypes.c_double(criteria_Dm))
            
            number_of_doublets_strict = np.sum(doublet_or_not_strict, axis = 0)
            number_of_doublets_general = np.sum(doublet_or_not_general, axis = 0)
            output_s.append(number_of_doublets_strict*2/particle_numbers)
            output_g.append(number_of_doublets_general*2/particle_numbers)
            # one doublet has two RBCs so to calculate doublet fraction, we need multiply # of doublets by 2
            if (make_plot == 1):
                ax1.plot(np.array(list(range(timesteps - period)))*WriteProps, number_of_doublets_strict*2/particle_numbers, label = "{}Dm, {}t_rot, strict".format(criteria_Dm, criteria_T))
                ax1.plot(np.array(list(range(timesteps - period)))*WriteProps, number_of_doublets_general*2/particle_numbers, label = "{}Dm, {}t_rot, general".format(criteria_Dm, criteria_T))
    

    # Make plot
    ###########################################################################
    if (make_plot == 1):
        ax1.set_xlabel("timesteps", fontsize = 20)
        ax1.set_ylabel("doublet fraction", fontsize = 20)
        ax1.set_title("h=24, phi={}, Re=0.1, Ca={}, D={}, eqWCA={}".format(phi, Ca, D, eqWCA), fontsize = 24)
        ax1.legend(fontsize = 16)
        ax2 = ax1.twinx()
        mn, mx = ax1.get_ylim()
        ax2.set_ylim(mn*particle_numbers/2, mx*particle_numbers/2)
        ax2.set_ylabel('# of doublets', fontsize = 20)
        #plt.savefig("./Pictures/Doublet_Fraction_{}_multiple_criteria.png".format(job_name), dpi = 300)
        #plt.close()
        plt.show()

    return [output_s, output_g, timesteps]


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
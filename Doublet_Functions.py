import numpy as np 
import matplotlib.pyplot as plt
import sys
import time
from scipy import fft, signal, stats

def calcDoubletFraction(input_phi, input_Ca, input_criteria_Ts, input_criteria_Dms, make_plot, *input_others):
    """
    Input:
    input_criteria_Ts and input_criteria_Dms must be either list or a float
    make_plot = 1 means making plots, 0 means doesn't need to make plots
    
    Output:
    A list of which element is the dounelt fraction time series for each set of criteria (T, r)
    """

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
    # Read the parameters
    with open("/Users/andrewliu/remote_disk/Data_Transfer/{}_parameter.txt".format(job_name)) as f:
        pre_parameters = f.readlines()
    timesteps = int((pre_parameters[pre_parameters.index("timesteps\n")+1])[:-1])
    particle_numbers = int((pre_parameters[pre_parameters.index("particle_numbers\n")+1])[:-1])
    interval = int((pre_parameters[pre_parameters.index("interval\n")+1])[:-1])
    points_per_particle = int((pre_parameters[pre_parameters.index("points_per_particle\n")+1])[:-1])
    
    # Two-cell system
    try:
        dim = eval((pre_parameters[pre_parameters.index("dim\n")+1])[:-1])
    # Suspension system
    except:
        dim = [144, 24, 144]

    # Read the position of center of mass
    with open("/Users/andrewliu/remote_disk/Data_Transfer/{}_COMs.txt".format(job_name)) as f:
        pre_COMs = eval(f.readlines()[0])
    COMs = np.reshape(pre_COMs, (particle_numbers, timesteps, 3))    

    # Read the y positions
    with open("/Users/andrewliu/remote_disk/Data_Transfer/{}_Ypos_t.txt".format(job_name)) as f:
        pre_Ypos_t = eval(f.readlines()[0])
    Ypos_t = np.reshape(pre_Ypos_t, (interval, particle_numbers*points_per_particle))
    

    
    # Calculate t_rot
    ###########################################################################
    # Format of Ypos_t is Ypos_t[t, particle_id]

    Periods = np.zeros(particle_numbers*points_per_particle)
    for i in range(particle_numbers*points_per_particle):
        Ypos_t_norm = Ypos_t[:, i] - np.mean(Ypos_t[:, i])
        f, Pxx = signal.periodogram(Ypos_t_norm, fs = 1, window='hanning', scaling='spectrum')
        valid_range = np.where(f > 0.01)
        Periods[i] = 1/(f[valid_range][np.argsort(Pxx[valid_range])[-1]])
    
    # Abandon arithmetic mean because there might be some sigular periods calculated by fft   
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
    diffpos = np.zeros((number_of_pairs, timesteps))
    indice_pairs = np.zeros((number_of_pairs, 2))

    count = 0
    for i in range(particle_numbers-1):
        for j in range(i+1, particle_numbers):
            indice_pairs[count, 0], indice_pairs[count, 1] = i, j
            diffpos[count, :] = np.linalg.norm((COMs[i, :, :] - COMs[j, :, :]), axis=1)
            count += 1       


    # Correct diffpos here
    for k in range(number_of_pairs):
        for t in range(1, timesteps):
            # exceeds maximum physical displacement, i.e. one of the two RBCs cross the boundry
            if abs(diffpos[k, t-1] - diffpos[k, t]) > 2*Dm: 
                i, j = int(indice_pairs[k,0]), int(indice_pairs[k,1])
                correct_current_pos = COMs[j, t, :] # modify the position of the latter one of the two RBCs

                # modify the x coordinate 
                if (COMs[i, t, 0]-int(dim[0]/2))*(COMs[j, t, 0]-int(dim[0]/2)) < 0:
                    if (COMs[j, t, 0]-int(dim[0]/2)) < 0: correct_current_pos[0] += int(dim[0])
                    else: correct_current_pos[0] -= int(dim[0])

                # modify the z coordinate 
                if (COMs[i, t, 2]-int(dim[1]/2))*(COMs[j, t, 2]-int(dim[1]/2)) < 0:
                    if (COMs[j, t, 2]-int(dim[1]/2)) < 0: correct_current_pos[2] += int(dim[1])
                    else: correct_current_pos[2] -= int(dim[1])

                # calculate the correct diff COM distance here
                diffpos[k, t] = np.linalg.norm(COMs[i, t, :] - correct_current_pos)


                
    # Calculate doublet fraction
    ###########################################################################

    if (make_plot == 1):
        fig, ax1 = plt.subplots(figsize = (16, 12))

    output = []
    for criteria_Dm in criteria_Dms:
        if criteria_Dm == 0.75: marker_Dm = '--'
        elif criteria_Dm == 1.0: marker_Dm = '-'

        for criteria_T in criteria_Ts:
            period = int(round(criteria_T*rotation_time))
            doublet_or_not = np.zeros((number_of_pairs, timesteps - period)) # 1 means there is a doublet then
            for i in range(number_of_pairs):# use KMP like algorithm to accelerate the calculation
                # initialize here
                t = 0
                try: # doesn't form doublet
                    t += ((np.argwhere(diffpos[i,0:period] > (criteria_Dm*Dm)).max())+1)
                except: # form doublet
                    doublet_or_not[i, 0] = 1
                    t += 1 # t = 1
                while t < timesteps - period:
                    if doublet_or_not[i, t-1] == 1: # this pair forms doublet in previous timestep
                        if diffpos[i, t+period-1] < criteria_Dm*Dm: # the doublet survives
                            doublet_or_not[i, t] = 1
                            t += 1
                        else: # the doublet breaks
                            t += period
                    else: # this pair doesn't form doublet in previous timestep
                        try: # doesn't form doublet
                            t += ((np.argwhere(diffpos[i,t:t+period] > (criteria_Dm*Dm)).max())+1)
                        except: # form doublet
                            doublet_or_not[i, t] = 1
                            t += 1 # t = 1
            number_of_doublets = np.sum(doublet_or_not, axis = 0)
            output.append(number_of_doublets*2/particle_numbers)
            # one doublet has two RBCs so to calculate doublet fraction, we need multiply # of doublets by 2
            if (make_plot == 1):
                ax1.plot(np.array(list(range(timesteps - period)))*WriteProps, number_of_doublets*2/particle_numbers, 
                        label = "{}Dm, {}t_rot".format(criteria_Dm, criteria_T), linestyle = marker_Dm)

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
        plt.close()

    return output


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
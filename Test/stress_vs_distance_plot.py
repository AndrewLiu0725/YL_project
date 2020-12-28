import numpy as np 
import matplotlib.pyplot as plt
from Doublet_Functions import calcDoubletFraction, getInstrinsicViscosity, getStress
import scipy as scipy

"""
This code is to plot elastic stress tensor vs interparticle distance, avg or max.
Also provided linear regression line.
"""

def make_patch_spines_invisible(ax):
    ax.set_frame_on(True)
    ax.patch.set_visible(False)
    for sp in ax.spines.values():
        sp.set_visible(False)

def plot():
    result = calcDoubletFraction(phi, Ca, 1.0, r, 0, angle)
    dfg = result[0][0]
    diffpos = result[3][0,0:len(dfg)]
    stress = getStress(phi, Ca, 0, 0, ncycle, angle)[:,3][0:len(dfg)]
    time = np.array(list(range(len(dfg))))

    # Consider phase shift
    w = 3
    st, et = 0, 1900
    local_min_time_distance = []
    local_min_time_stress = []

    avg_distance = []
    max_distance = []
    avg_stress = []
    max_stress = []

    distance_sum, distance_period, distance_max = 0, 0, -1
    stress_sum, stress_period, stress_max = 0, 0, -1

    min_threshold = 2 # for 1 cycle
    num_min_distance, distance_min_num = 0, 0    

    prev_dis_min_time = 0
    # Find local min time for d_sep
    for t in time[st+w:et-w]:
        # in doubelt state
        if dfg[t] == 1: 
            # distance
            # local min
            if ((((diffpos[t-1]-diffpos[t]) > 0) and ((diffpos[t+1]-diffpos[t]) > 0)) and
                (((diffpos[t-3]-diffpos[t]) > 0) and ((diffpos[t+3]-diffpos[t]) > 0))):
                if num_min_distance < min_threshold:
                    num_min_distance += 1
                else:
                    # Mark the starting and ending time of this cycle
                    local_min_time_distance.append(prev_dis_min_time)
                    local_min_time_distance.append(t)
                    distance_min_num += 1
                    avg_distance.append(distance_sum/distance_period)
                    max_distance.append(distance_max)
                distance_sum, distance_period, distance_max = 0, 0, -1
                prev_dis_min_time = t

            # interval
            else:
                distance_period += 1
                distance_sum += diffpos[t]
                # max
                if diffpos[t] > distance_max:
                    distance_max = diffpos[t]
        else:
            num_min_distance = 0


    # Find the corresponding local min time for stress tensor
    st_or_et = 0 # 0 means finding st, 1 means finding et
    stress_min_num = 0
    prev_stress_min_time = 0
    for t in time[st+w:et-w]:
        # stress
        # local min
        if ((((stress[t-1]-stress[t])>0) and ((stress[t+1]-stress[t])>0)) and
            (((stress[t-3]-stress[t])>0) and ((stress[t+3]-stress[t])>0))):
            # True min
            # find st
            if st_or_et == 0:
                if t >= local_min_time_distance[2*stress_min_num]:
                    # cur min is st
                    if (t <= (local_min_time_distance[2*stress_min_num]+1)): # +2 for error tolerance
                        if (stress[t-1] - stress[t]) > 0.0001*0.8:
                            local_min_time_stress.append(prev_stress_min_time)
                        else:
                            local_min_time_stress.append(t)
                    # prev min is st
                    else:
                        local_min_time_stress.append(prev_stress_min_time)
                    st_or_et = 1
            # find et
            if st_or_et == 1:
                if (t-local_min_time_stress[2*stress_min_num]) > 0.6*(local_min_time_distance[2*stress_min_num+1]-local_min_time_distance[2*stress_min_num]):
                    local_min_time_stress.append(t)
                    stress_min_num += 1
                    if stress_min_num == distance_min_num: break
                    st_or_et = 0
                    # set the next st if possible
                    if local_min_time_distance[2*stress_min_num] == local_min_time_distance[2*stress_min_num-1]:
                        local_min_time_stress.append(t)
                        st_or_et = 1
            prev_stress_min_time = t
                    
    # Calc the avg and max of stress tensor according to the local time array created in the previous step
    for i in range(stress_min_num):
        avg_stress.append(np.mean(stress[local_min_time_stress[2*i] : local_min_time_stress[2*i+1]]))
        max_stress.append(np.max(stress[local_min_time_stress[2*i] : local_min_time_stress[2*i+1]]))
    
    
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))
    fig.suptitle('phi = {}, Ca = {}, angle = {}, criteria r = {}Dm\n'.format(phi, Ca, angle, r), fontsize = 15)
    slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(avg_distance, avg_stress)
    fit_dis = np.array(avg_distance)
    fit_stress = fit_dis*slope + intercept
    

    ax1.set_title("Consider phase shift")
    ax1.plot(avg_distance, avg_stress, 'o', label = "Data")
    ax1.plot(fit_dis, fit_stress, label = "Fitted line")
    ax1.legend()
    ax1.annotate("y={}x+{}\n".format(np.round(slope,7), np.round(intercept,6)) + r"$R^{2}$" + "={}".format(np.round(r_value**2,6)), xy=(0.55, 0.75), xycoords='axes fraction')
    ax1.set_xlabel(r'$d _{sep}$', fontsize = 15)
    ax1.set_ylabel(r'$\sigma _{elas, yx}$', fontsize = 15)

    # Not consider phase shift
    w = 3
    st, et = 0, 1900

    avg_distance = []
    max_distance = []
    avg_stress = []
    max_stress = []

    num_min = 0
    threshold_distance = 0.5
    threshold_stress = 0.0001
    distance_sum, distance_period, distance_max = 0, 0, -1
    stress_sum, stress_max = 0, -1

    min_threshold = 1
    num_min_stress, stress_min_num = 0, 0

    for t in time[st+w:et-w]:
        if dfg[t] == 1: # in doubelt state
            # interval
            if num_min >= min_threshold: # start storing values
                distance_period += 1
                # avg
                distance_sum += diffpos[t]
                stress_sum += stress[t]
                # max
                if diffpos[t] > distance_max:
                    distance_max = diffpos[t]
                if stress[t] > stress_max:
                    stress_max = stress[t]


            # local min
            if (((diffpos[t-1]-diffpos[t]) > 0) and ((diffpos[t+1]-diffpos[t]) > 0)):
                if (((diffpos[t-w]-diffpos[t]) > threshold_distance) and ((diffpos[t+w]-diffpos[t]) > threshold_distance)):
                    if num_min < min_threshold:
                        num_min += 1
                    else:
                        avg_distance.append(distance_sum/distance_period)
                        max_distance.append(distance_max)
                        avg_stress.append(stress_sum/distance_period)
                        max_stress.append(stress_max)
        
                        distance_sum, distance_period, distance_max = 0, 0, -1
                        stress_sum, stress_max = 0, -1

            if (((stress[t-1]-stress[t])>0) and ((stress[t+1]-stress[t])>0)):
                if (((stress[t-w]-stress[t])>threshold_stress) and ((stress[t+w]-stress[t])>threshold_stress)):
                    if num_min_stress < min_threshold:
                        num_min_stress += 1
                    else:
                        stress_min_num += 1

        else:
            num_min = 0
            num_min_stress = 0

    slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(avg_distance, avg_stress)
    fit_dis = np.array(avg_distance)
    fit_stress = fit_dis*slope + intercept

    ax2.set_title("Not consider phase shift")
    ax2.plot(avg_distance, avg_stress, 'o', label = "Data")
    ax2.plot(fit_dis, fit_stress, label = "Fitted line")
    ax2.legend()
    ax2.annotate("y={}x+{}\n".format(np.round(slope,7), np.round(intercept,6)) + r"$R^{2}$" + "={}".format(np.round(r_value**2,6)), xy=(0.55, 0.75), xycoords='axes fraction')
    ax2.set_xlabel(r'$d _{sep}$', fontsize = 15)
    ax2.set_ylabel(r'$\sigma _{elas, yx}$', fontsize = 15)
    #plt.tight_layout()
    plt.show()
    
    '''
    if flag == 2:
        fig, host = plt.subplots(figsize = (12, 9))
        fig.subplots_adjust(right=0.75)

        par1 = host.twinx()
        par2 = host.twinx()
        par2.spines["right"].set_position(("axes", 1.1))
        make_patch_spines_invisible(par2)
        par2.spines["right"].set_visible(True)

        simul_time = list(range(len(dfg)))
        p1, = host.plot(simul_time, dfg, "b-", label = "doublet fraction")
        p2, = par1.plot(simul_time, diffpos, "r-", label = "interparticle distance")
        par1.plot(local_min_time_distance, diffpos[local_min_time_distance], "x")
        p3, = par2.plot(simul_time, stress, "g-", label = r'$\sigma _{elas, yx}$')
        par2.plot(local_min_time_stress, stress[local_min_time_stress], "o")

        host.set_title('phi = {}, Ca = {}, angle = {}\ncriteria r = {}Dm'.format(phi, Ca, angle, r), fontsize = 30)
        host.set_xlabel("timesteps(strain)", fontsize = 20)
        host.set_ylabel("Doublet fraction", fontsize = 20)
        host.set_ylim([-0.1, 1.1])
        host.set_xlim([le, re])
        par1.set_ylabel("interparticle distance", fontsize = 20)
        par2.set_ylabel(r'$\sigma _{elas, yx}$', fontsize = 20)

        host.yaxis.label.set_color(p1.get_color())
        par1.yaxis.label.set_color(p2.get_color())
        par2.yaxis.label.set_color(p3.get_color())

        tkw = dict(size=6, width=1.5)
        host.tick_params(axis='y', colors=p1.get_color(), **tkw)
        par1.tick_params(axis='y', colors=p2.get_color(), **tkw)
        par2.tick_params(axis='y', colors=p3.get_color(), **tkw)
        host.tick_params(axis='x', **tkw)

        lines = [p1, p2, p3]
        host.legend(lines, [l.get_label() for l in lines], prop={'size': 12}, bbox_to_anchor=(1.2,1.2))
        fig.tight_layout()
        plt.show()
    '''

ncycle = 2000
r = 1
phi, Ca, angle = 6.9, 0.05, -30
le, re = 200, 300
plot()
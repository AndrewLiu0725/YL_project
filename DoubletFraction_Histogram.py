# ===============================================================================
# Copyright 2021 An-Jun Liu
# Last Modified Date: 01/30/2021
# ===============================================================================
import numpy as np 
import matplotlib.pyplot as plt
import sys
import time
from RBC_Utilities import calcDoubletFraction
import os 
from scipy import stats

"""
This code is to plot doublet fraction's histogram and fit it with skew normal distribution in suspension system
"""

start_time = time.time()

path = "/raid6/ctliao/Data/HI_ordering/"

def rmse(predictions, targets):
    return np.sqrt(((predictions - targets) ** 2).mean())

# parser
def parser(string):
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



r = 1

def makeDFHist(result):
    for i in range(2): # for whole and second half time series
        if i == 0:
            data = (result[0][0]*result[4]/2).astype(int) # convert to number of doublets
        else:
            data = (result[0][0][-int(result[1]/2):]*result[4]/2).astype(int) # convert to number of doublets
        
        fig = plt.figure(figsize=(8, 6))
        ax1 = fig.add_subplot(111)    
        ax2 = ax1.twiny()
        max_num_of_doublet = max(data)
        min_num_of_doublet = min(data)
        val_original, cnt_original = np.unique(data, return_counts=True) # the frequency of each value

        # fill with zeros
        val, cnt = np.zeros(max_num_of_doublet-min_num_of_doublet+1), np.zeros(max_num_of_doublet-min_num_of_doublet+1)
        ptr = 0
        for index, numDoublet in enumerate(val_original):
            val[ptr] = numDoublet
            cnt[ptr] = cnt_original[index]
            ptr += 1
            if index == (len(val_original)-1): break # end
            if (val_original[index+1] > numDoublet + 1):
                for i in range(numDoublet + 1, val_original[index+1]):
                    val[ptr] = i
                    cnt[ptr] = 0 # fill with zero to calcuate RMSD
                    ptr += 1
                    
        pmf = cnt / len(data) # convert to probability mass function
        ax1.bar(val, pmf, label = "Simulation data") # plot pmf

        if (len(pmf.nonzero()[0]) == 1):
            pass
        else:
            # fit with skew normal
            a, loc, scale = stats.skewnorm.fit(data)
            x = np.linspace(min_num_of_doublet, max_num_of_doublet, 100)
            p = stats.skewnorm.pdf(x, a, loc, scale)
            ax1.plot(x, p, 'k', linewidth=2, label = "Fitted line")
            # calculate rms error
            x_predicted = np.array(list(range(min_num_of_doublet, max_num_of_doublet + 1)))
            df_estimated = stats.skewnorm.pdf(x_predicted, a, loc, scale)
            rms = rmse(df_estimated, pmf)
            ax1.text(max_num_of_doublet*0.85, max(pmf)*0.6, "RMSD = {:.3}\na = {:.3}\nloc = {:.3}\nscale = {:.3}".format(rms, a, loc, scale), weight = "bold")

        # set second x axes
        def tick_function(x):
            return 2*x/result[4]
        ax2.set_xbound(tick_function(ax1.get_xbound()[0]), tick_function(ax1.get_xbound()[1]))
        ax1.set_xlabel("# of doublets", fontsize = 15) 
        ax1.set_ylabel("Frequency", fontsize = 15)
        if i == 0:
            ax1.set_title("Doublet fraction's histogram\nphi = {}, Ca = {} (whole time series)".format(phi, Ca), fontsize = 20)
        else:    
            ax1.set_title("Doublet fraction's histogram\nphi = {}, Ca = {} (second half time series)".format(phi, Ca), fontsize = 20)
        ax1.legend()
        ax2.set_xlabel('Doublet fraction', fontsize = 15)
        fig.tight_layout()
        job_name = "h24phi{}Re0.1Ca{}WCA1zero0.8-{}".format(phi, Ca, ensemble_id)
        if i == 0:
            plt.savefig("./Pictures/hist/{}_histogram_wholetimeseries.png".format(job_name), dpi = 200)
        else:
            plt.savefig("./Pictures/hist/{}_histogram_secondhalftimeseries.png".format(job_name), dpi = 200)
        plt.close()


# run over all sets of parameters
for phi in parameter_set.keys():
    for Ca in parameter_set[phi].keys():
        for ensemble_id in parameter_set[phi][Ca]:
            try:
                result = calcDoubletFraction(phi, Ca, 1, r, 0, ensemble_id, 1)
                makeDFHist(result)
            except:
                pass
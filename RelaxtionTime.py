# ===============================================================================
# Copyright 2021 An-Jun Liu
# Last Modified Date: 03/20/2021
# ===============================================================================
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import time
import datetime
import pickle
import random
from statsmodels.tsa.stattools import adfuller

"""
This code is to calculate the relaxation time for twoc-ell and suspension system
"""

DEBUG = 0

start_time = time.time()

C_0 = 0 # global C(0)

def C_tau_twoVar(x, a, b):
    # a = C(0), b = t_relax
    return a*np.exp(-x/b)

def C_tau_oneVar(x, a):
    # a = t_relax
    #if DEBUG: print("In C(tau), C_0 = ", C_0)
    return C_0*np.exp(-x/a)

def calcRelaxationTime(phi, Ca, significance_level, plot):
    """
    Input:

    significance_level: alpha value in null hypothesis test, recommeded value: 0.05

    func_type: 1 is for two-var and 0 is for one-var C(tau)


    Output:

    [Is steady state?, t_relax, standard deviation errors on the parameters]
    """

    df = avg_df_dict[phi][Ca]

    # calculate stationarity
    alpha = significance_level
    starting_point = [0, 0.25, 0.5]
    p_value = np.zeros(3)
    unstable = 0

    for i in range(3):
        result = adfuller(df[int(starting_point[i]*len(df)):])
        p_value[i] = result[1]

        # p value > alpha, fail to reject the null hypothesis (non-stationary)
        # means the time series is not stationary
        if p_value[i] > alpha:
            unstable = 1

    if DEBUG: print(p_value)


    # calculate the relxation time

    # setup
    steady_f = np.mean(df[-int(0.25*len(df)):]) # the avg of the last quarter time series
    deviation_f = df - steady_f
    T = len(df)
    if DEBUG: print("T = ", T)
    tau_range = np.array([i for i in range(0, T)])
    C = []

    # calculate C(tau)
    for tau in tau_range:
        C_tmp = np.dot(deviation_f[:T-tau], deviation_f[tau:])/(T-tau)
        C.append(C_tmp)

    # fitting
    t_relax, rmse = [], []
    C = np.array(C)
    global C_0
    C_0 = C[0]
    
    end_list = [T, int(T/2), int(T/4)]
    text_list = ["whole", "half", "quarter"]
    color_list = ["r-", "g-", "b-"]

    if plot:
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize = (9, 6))
        fig.suptitle("Relaxation time ({}): phi = {}, Ca = {}".format(system, phi, Ca), fontsize = 20)
        ax1.plot(list(range(T)), df)
        ax1.set_title("Doublet fraction time series", fontsize = 15)
        ax1.set(xlabel = "time", ylabel = "doublet fraction")

        ax2.plot(tau_range, C, label = "original data")
        ax2.set_title(r'$C\left( \tau \right)$')
        ax2.set(xlabel = r'$\tau$', ylabel = "correlation")


    for i in range(3):
        popt, pcov = curve_fit(C_tau_oneVar, tau_range[:end_list[i]], C[:end_list[i]])
        #perr = np.sqrt(np.diag(pcov))

        if plot:
            ax2.plot(tau_range, C_tau_oneVar(tau_range, *popt), color_list[i], label = text_list[i]+"\nt_relax={}".format(int(popt[0])))

        t_relax.append(popt[0])
        rmse.append(np.sqrt(np.mean(np.square(C-(C_tau_oneVar(tau_range, *popt))))))


    if plot:
        ax2.legend()
        fig.tight_layout()
        plt.subplots_adjust(top = 0.85)
        plt.show()

    return [False, None, None] if unstable else [True, t_relax, rmse]

system = "TwoCell"

# suspension part
with open("Data/AverageDF_{}.pickle".format(system), 'rb') as handle:
    avg_df_dict = pickle.load(handle)

TEST = 0

if TEST:
    #print(calcRelaxationTime(avg_df_dict[4.9488][0.08], 0.05, 500, 1, 1))
    #print(calcRelaxationTime(4.3, 0.13, 0.05, 1))
    #print(calcRelaxationTime(4.0, 0.1, 0.05, 1))
    
    
    for _ in range(5):
        phi = random.choice(list(avg_df_dict.keys()))
        Ca = random.choice(list(avg_df_dict[phi].keys()))
        print(calcRelaxationTime(phi, Ca, 0.05, 1))
    
    

else:
    relaxation_time = []
    total_count = 0
    for phi in avg_df_dict.keys():
        for Ca in avg_df_dict[phi].keys():
            total_count += 1
            result = calcRelaxationTime(phi, Ca, 0.05, 0)
            if result[0]:
                relaxation_time.append(np.mean(result[1]))

    plt.hist(relaxation_time)
    plt.xlabel("Relaxation time")
    plt.ylabel("Count")
    plt.title("Relaxation time ({}, {}%)".format(system, round(len(relaxation_time)/total_count, 2)*100))
    #plt.show()
    plt.savefig("./Pictures/{}System_RelaxationTimeHistogram.png".format(system), dpi = 200)

print('\nTotal time elapsed = {}'.format(str(datetime.timedelta(seconds=time.time()-start_time))))
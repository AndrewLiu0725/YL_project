# ===============================================================================
# Copyright 2021 An-Jun Liu
# Last Modified Date: 05/12/2021
# ===============================================================================
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import time
import datetime
import pickle
from statsmodels.tsa.stattools import adfuller
import random

"""
This code is to calculate the relaxation time for two-cell and suspension system
"""

DEBUG = 1

C_0 = 0 # global C(0)

def C_tau_oneVar(x, a):
    # a = t_relax
    return C_0*np.exp(-x/a)

def calcRelaxationTime(ts, significance_level, maxlag_p, plot):
    """
    Input:

    ts: 1-D numpy array

    significance_level: alpha value in null hypothesis test, recommeded value: 0.05

    maxlag_p: maximum lag which is included in test. 
    The adfuller function will automatically determining the lag length among the values [0, maxlag].
    The default method is that the number of lags is chosen to minimize the corresponding information criterion

    Output:

    [Is steady state?, t_relax, standard deviation errors on the parameters]
    """
    
    ### calculate stationarity
    T = len(ts)
    unstable = 0
    fail = 1

    try:
        adf_result = adfuller(ts[int(T/2):], maxlag = maxlag_p)
        if DEBUG: print("p-value = ",adf_result[1])
        # p value > alpha, fail to reject the null hypothesis (non-stationary) means the time series is not stationary
        if adf_result[1] >= significance_level:
            unstable = 1
        # p value <= alpha, the last test_length time series is stationary according to the adf test
        else:
            fail = 0
    
    except Exception as e:
        print(e)


    ### calculate the relxation time
    # setup
    steady_ts = np.mean(ts[int(T/2):]) # the avg of the test_length time series
    deviation_ts = ts - steady_ts
    tau_range = np.arange(0, T)
    C = np.zeros(T)

    # calculate C(tau)
    for tau in tau_range:
        C[tau] = np.dot(deviation_ts[:T-tau], deviation_ts[tau:])/(T-tau)

    # fitting
    t_relax, rmse = [], []
    global C_0 # mark that the C_0 is a global variable (must)
    C_0 = C[0]
    
    end_list = [T, int(T/2), int(T/4)]
    text_list = ["whole", "half", "quarter"]
    color_list = ["r-", "g-", "b-"]

    slicing = 10
    if plot:
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize = (9, 6))
        fig.suptitle("Relaxation Time Analysis", fontsize = 20)
        ax1.plot(tau_range[::slicing], ts[::slicing])
        ax1.set_title("Time series", fontsize = 15)
        ax1.set(xlabel = "time", ylabel = "variable")

        ax2.plot(tau_range, C, label = "original data")
        ax2.set_title(r'$C\left( \tau \right)$')
        ax2.set(xlabel = r'$\tau$', ylabel = "correlation")


    for i in range(3):
        try:
            popt = curve_fit(C_tau_oneVar, tau_range[:end_list[i]], C[:end_list[i]])[0]

            if plot:
                ax2.plot(tau_range, C_tau_oneVar(tau_range, *popt), color_list[i], label = text_list[i]+"\nt_relax={}".format(int(popt[0])))

            t_relax.append(popt[0])
            rmse.append(np.sqrt(np.mean(np.square(C-(C_tau_oneVar(tau_range, *popt))))))
        except:
            pass


    if plot:
        ax2.legend()
        fig.tight_layout()
        plt.subplots_adjust(top = 0.85)
        plt.show()

    return [False] if (unstable or fail) else [True, t_relax, rmse]


def analyzeRelaxationTime(variable_type, system, time_series_type, significance_level, plot, save):
    """
    Input:

    plot: flag to make histogram of relaxation time for stationary cases

    save: flag to store unstationary ratio
    """
    # get the nested dictionary of time series
    with open("{}_{}_{}.pickle".format(variable_type, system, time_series_type), 'rb') as handle:
        ts_dict = pickle.load(handle)

    # collect the relaxation time
    relaxation_time = []
    total_count = 0
    autoregressive_order_p = 20
    unstationary_ratio = []

    for phi in ts_dict.keys():
        for Ca in ts_dict[phi].keys():
            unstationary_count = 0
            ensemble_count = 0

            # run over ensemble
            for ts in ts_dict[phi][Ca]:
                total_count += 1
                ensemble_count += 1
                result = calcRelaxationTime(ts, significance_level, autoregressive_order_p, 0)
                if result[0]:
                    relaxation_time.append(np.mean(result[1]))
                else:
                    unstationary_count += 1
            
            unstationary_ratio.append([phi, Ca, unstationary_count/ensemble_count])

    if plot:
        plt.hist(relaxation_time)
        plt.xlabel("Relaxation time ({})".format(r'$\dot \gamma t$'))
        plt.ylabel("Count")
        plt.title("Relaxation Time Histogram\n({}, {}, {}\nsignificance level = {}, maxlag = {}, {}%)".format("Intrinsic Viscosity" if variable_type == "IV" else "Doublet Fraction", 
        system, time_series_type, significance_level, autoregressive_order_p, round(len(relaxation_time)/total_count, 2)*100))
        plt.savefig("./Pictures/{}System_{}_RelaxationTime_{}_alpha_{}_p_{}_Histogram.png".format(system, "IntrinsicViscosity" if variable_type == "IV" else "DoubletFraction", 
        time_series_type, significance_level, autoregressive_order_p), dpi = 200)
        plt.close()

    if save:
        with open("UnstationaryRatio_{}_{}_{}_alpha_{}_p_{}.pickle".format(variable_type, system, time_series_type, significance_level, autoregressive_order_p), 'wb') as handle:
            pickle.dump(unstationary_ratio, handle, protocol=pickle.HIGHEST_PROTOCOL)


def test(variable_type, system, time_series_type, random_flag, phi, Ca):   
    """
    random flag: 1 for randomly choosing [phi, Ca], 
    """
    with open("./Data/{}_{}_{}.pickle".format(variable_type, system, time_series_type), 'rb') as handle:
        ts_dict = pickle.load(handle)
    autoregressive_order_p = 50
    if random_flag:
        phi = random.choice(list(ts_dict.keys()))
        Ca = random.choice(list(ts_dict[phi].keys()))
    print("phi = {}, Ca = {}".format(phi, Ca))
    print(calcRelaxationTime(ts_dict[phi][Ca][0][1:], 0.05, autoregressive_order_p, 1000, 1))


if __name__ == "__main__":
    start_time = time.time()

    system_list = ["TwoCell", "Suspension"]
    #time_series_type_list = ["Indivisual", "EnsembleAveraged"]
    time_series_type_list = ["EnsembleAveraged"]

    for system in system_list:
        for time_series_type in time_series_type_list:
            for alpha in [0.05, 0.1]:
                print("Current task: {}, {}, alpha = {}\n".format(system, time_series_type, alpha))
                analyzeRelaxationTime("IV", system, time_series_type, alpha, 1, 1)

    print('\nTotal time elapsed = {}'.format(str(datetime.timedelta(seconds=time.time()-start_time))))
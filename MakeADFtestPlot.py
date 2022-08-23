# ===============================================================================
# Copyright 2021 An-Jun Liu
# Last Modified Date: 08/22/2022
# ===============================================================================
import pickle
import matplotlib.pyplot as plt 
from matplotlib.ticker import MaxNLocator
import numpy as np 
from statsmodels.tsa.stattools import adfuller

plt.rcParams['font.family'] = 'DeJavu Serif'
plt.rcParams['font.serif'] = ['Times New Roman']

def makePlot(phi, Ca, system, ax:plt.Axes):
    with open("./Data/IV_{}_EnsembleAveraged.pickle".format(system), 'rb') as handle:
        ts_dict = pickle.load(handle)
    
    ts = ts_dict[phi][Ca][0]
    significance_level = 0.1

    ### calculate stationarity
    T = len(ts)
    unstable = 0
    fail = 1
    adf_result = adfuller(ts[int(T/2):], maxlag = 20)
    # p value > alpha, fail to reject the null hypothesis (non-stationary) means the time series is not stationary
    if adf_result[1] >= significance_level:
        unstable = 1
    # p value <= alpha, the last test_length time series is stationary according to the adf test
    else:
        fail = 0

    slicing = 1
    factor = 1 if system == "Two-Cell" else 4000/3669

    ax.plot(np.arange(T)[::slicing], ts[::slicing])
    ax.set_xlabel("{}".format(r'$\gamma _{rot}$'), fontsize = 20)
    ax.set_ylabel(r'$\left[ \eta \right]$', fontsize = 20)
    ax.set_title("{} System: {} = {}, Ca = {}\np-value = {}".format(system, r'$\phi$', phi, Ca, adf_result[1]), fontsize = 15)
    ax.tick_params(which='both', labelsize=15, width=2, length=8, direction='in')
    ax.xaxis.set_major_locator(MaxNLocator(5))
    ax.xaxis.set_minor_locator(MaxNLocator(10))
    ax.yaxis.set_major_locator(MaxNLocator(5)) 
    #ax.yaxis.set_minor_locator(MaxNLocator(10))


# main code
# ===============================================================================
rows, columns = 2, 2
fig, axs = plt.subplots(rows, columns, figsize = (7.5*columns, 6*rows))

parameter_list = [[2.9993, 0.14, 'Suspension'], [3.8991, 0.06, 'Suspension'], [4.0, 0.01, 'TwoCell'], [4.7, 0.07, 'TwoCell']]
for i in range(4):
    makePlot(parameter_list[i][0], parameter_list[i][1], parameter_list[i][2], axs[i//2, i%2])

fig.tight_layout()
fig.savefig("Pictures/ADFtest/fig_S6.png", dpi = 200)
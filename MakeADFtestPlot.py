# ===============================================================================
# Copyright 2021 An-Jun Liu
# Last Modified Date: 06/10/2021
# ===============================================================================
import pickle
with open("./Data/IV_Suspension_EnsembleAveraged.pickle", 'rb') as handle:
    ts_dict = pickle.load(handle)

import matplotlib.pyplot as plt 
import numpy as np 
from statsmodels.tsa.stattools import adfuller

phi = 3.8991
Ca = 0.06
system = "Suspension"
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
fig = plt.figure(figsize = (7.5, 6))
plt.plot(np.arange(T)[::slicing], ts[::slicing])
plt.title("Intrinsic Viscosity Time Series\n{} System: {} = {}, Ca = {}\np-value = {}".format(system, r'$\phi$', phi, Ca, adf_result[1]), fontsize = 15)
plt.xlabel("time ({})".format(r'$\dot \gamma t$'), fontsize = 15)
plt.ylabel(r'$\eta _{int}$', fontsize = 15)
#plt.show()
plt.tight_layout()
plt.savefig("Pictures/ADFtest/{}System_IV_phi_{}_Ca_{}.png".format("TwoCell" if system == "Two-Cell" else "Suspension", phi, Ca), dpi = 200)
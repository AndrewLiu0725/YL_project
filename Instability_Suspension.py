# ===============================================================================
# Copyright 2021 An-Jun Liu
# Last Modified Date: 03/13/2021
# ===============================================================================
import numpy as np
import matplotlib.pyplot as plt
import time
import datetime
from RBC_Utilities import calcDoubletFraction, getIntrinsicViscosity, getSuspensionParameterSets
import os
import sys

"""
This code is to calculate the instability of each ensemble-averaged doublet fraction time series 
and make a histogram out of it for suspension system.
The instability is defined as (avg(fourth quarther) - avg(third quarter))/avg(third quarter).
"""

FLAG = 1 # 1 for intrinsic viscosity, 0 for doublet fraction
DEBUG = 0
PRINT = 0
HIST = 1

start_time = time.time()

[phis, parameter_set] = getSuspensionParameterSets()

# calculate deviation
r = 1

deviation = []
ps = []
max_timesteps = 10000

for phi in parameter_set.keys():
    if phi in phis:
        for Ca in parameter_set[phi].keys():
            # initialize
            ensemble_count = 0
            df = np.zeros((len(parameter_set[phi][Ca]), max_timesteps))
            min_timesteps = max_timesteps

            # run over ensemble
            for ensemble_id in parameter_set[phi][Ca]:
                try:
                    if FLAG:
                        iv = getIntrinsicViscosity(phi, Ca, max_timesteps, ensemble_id, 1)
                        min_timesteps = min(min_timesteps, len(iv))
                        df[ensemble_count, :len(iv)] = iv[:len(iv)]
                    else:
                        result = calcDoubletFraction(phi, Ca, 1, r, 0, ensemble_id, 1)
                        min_timesteps = min(min_timesteps, result[1])
                        df[ensemble_count, :result[1]] = result[0][0][:result[1]]
                    ensemble_count += 1

                except KeyboardInterrupt:
                    print("Pressed ctrl+C\n")
                    sys.exit()

                except Exception as e:
                    print("Error: Case (phi = {}, Ca = {}, ensemble_id = {}): \n{}".format(phi, Ca, ensemble_id, e))
                    #print("Error: no preprocessed data in (phi = {}, Ca = {}, ensemble id = {})".format(phi, Ca, ensemble_id))

            if DEBUG: print("(phi = {}, Ca = {}): min timesteps = {}".format(phi, Ca, min_timesteps))

            if (ensemble_count > 0):
                avg_df = np.mean(df[:ensemble_count, :min_timesteps], axis=0)
                first = np.mean(avg_df[int(0.5*min_timesteps):int(0.75*min_timesteps)])
                second = np.mean(avg_df[int(0.75*min_timesteps):min_timesteps])
                deviation.append(abs((second-first)/first))
                ps.append([phi, Ca])
            else:
                print("Error: no data in (phi = {}, Ca = {})".format(phi, Ca))
            

a = sorted(range(len(deviation)), key=lambda k: deviation[k], reverse = True)    

if PRINT:
    for i in a:
        print("({}, {}) deviation = {}".format(ps[i][0], ps[i][1], deviation[i]))   

if HIST:
    plt.hist(deviation, bins = 20, weights = (np.ones_like(deviation)/len(deviation)))
    plt.xlabel("Instability")
    plt.ylabel("Density")
    plt.title("Instability histogram (suspension system)")
    plt.savefig("./Pictures/SuspensionSystem_Instability_Histogram.png", dpi = 300)
    plt.close()

print('Total time elapsed = {}'.format(str(datetime.timedelta(seconds=time.time()-start_time))))
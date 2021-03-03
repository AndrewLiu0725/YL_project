# ===============================================================================
# Copyright 2021 An-Jun Liu
# Last Modified Date: 02/27/2021
# ===============================================================================
import numpy as np
import matplotlib.pyplot as plt
import time
import datetime
from RBC_Utilities import calcDoubletFraction, getIntrinsicViscosity
import os
import sys

"""
This code is to calculate the instability of each ensemble-averaged doublet fraction time series 
and make a histogram out of it for two-cell system.
The instability is defined as (avg(fourth quarther) - avg(third quarter))/avg(third quarter).
"""

FLAG = 1 # 1 for intrinsic viscosity, 0 for doublet fraction
DEBUG = 0
PRINT = 0
HIST = 1

start_time = time.time()

# Parameters
r = 1.0
ncycle = 2000

vol = 746.3163
phi_range =[]
for x in range(30, 41):
    phi = 2*vol/(24*x**2)
    phi_range.append(float(str(phi*100)[:3]))
Ca_range = [i*0.01 for i in range(1, 21)]
angle_range = [90 - 10*j for j in range(18)]

deviation_ensemble = []
deviation_all = []
ps_ensemble = []
max_timesteps = 10000 # used to calculate ensemble one
shift = 0

# calculate the instability
for phi in phi_range:
    for Ca in Ca_range:

        # initialize
        ensemble_count = 0
        df = np.zeros((len(angle_range), max_timesteps))
        min_timesteps = max_timesteps

        # run over ensemble
        for angle in angle_range:
            try:
                if FLAG:
                    iv = getIntrinsicViscosity(phi, Ca, ncycle, angle, 0)
                    min_timesteps = min(min_timesteps, len(iv))
                    df[ensemble_count, :len(iv)] = iv[:len(iv)]
                else:
                    result = calcDoubletFraction(phi, Ca, 1, r, 0, angle, 0)
                    min_timesteps = min(min_timesteps, result[1])
                    df[ensemble_count, :result[1]] = result[0][0][:result[1]]
                
                first = np.mean(df[ensemble_count, int(0.5*min_timesteps):int(0.75*min_timesteps)]) + shift
                second = np.mean(df[ensemble_count, int(0.75*min_timesteps):min_timesteps]) + shift
                deviation_all.append(abs((second-first)/first))

                ensemble_count += 1

            except KeyboardInterrupt:
                print("Pressed ctrl+C\n")
                sys.exit()

            except Exception as e:
                print("Error: Case (phi = {}, Ca = {}, angle = {}): \n{}".format(phi, Ca, angle, e))

        if DEBUG: print("(phi = {}, Ca = {}): min timesteps = {}".format(phi, Ca, min_timesteps))

        if (ensemble_count > 0):
            avg_df = np.mean(df[:ensemble_count, :min_timesteps], axis = 0)
            first = np.mean(avg_df[int(0.5*min_timesteps):int(0.75*min_timesteps)]) + shift
            second = np.mean(avg_df[int(0.75*min_timesteps):min_timesteps]) + shift
            deviation_ensemble.append(abs((second-first)/first))
            ps_ensemble.append([phi, Ca])
        else:
            print("Error: no data in (phi = {}, Ca = {})".format(phi, Ca))

a = sorted(range(len(deviation_ensemble)), key=lambda k: deviation_ensemble[k], reverse = True)    

if PRINT:
    for i in a:
        print("({}, {}) deviation = {}".format(ps_ensemble[i][0], ps_ensemble[i][1], deviation_ensemble[i]))     
if HIST:
    plt.hist(deviation_ensemble, bins = 20, weights = (np.ones_like(deviation_ensemble)/len(deviation_ensemble)))
    plt.xlabel("Instability")
    plt.ylabel("Density")
    plt.title("Instability histogram (two-cell system, angular ensembled)")
    plt.savefig("./Pictures/TwoCellSystem_Instability_Histogram_AngularEnsembled.png", dpi = 300)
    plt.close()

    plt.hist(deviation_all, bins = 20, weights = (np.ones_like(deviation_all)/len(deviation_all)))
    plt.xlabel("Instability")
    plt.ylabel("Density")
    plt.title("Instability histogram (two-cell system)")
    plt.savefig("./Pictures/TwoCellSystem_Instability_Histogram.png", dpi = 300)
    plt.close()

print('Total time elapsed = {}'.format(str(datetime.timedelta(seconds=time.time()-start_time))))
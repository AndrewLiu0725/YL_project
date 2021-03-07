# ===============================================================================
# Copyright 2021 An-Jun Liu
# Last Modified Date: 03/07/2021
# ===============================================================================
import os
import subprocess
import time
import datetime

'''
This program is a script for preprocessing the data for two-cell system
'''

start_time = time.time()

# parameters
nstrain = 4000
vol = 746.3163
test = 0 # 1 for creating init case for each set of parameters
Ca_list = [round(i*0.01+0.01, 5) for i in range(20)]
phi_list = [float(str(2*vol/(24*x**2)*100)[:3]) for x in range(30, 41)] # y is fixed to 24
angle_list = [90 - 10*i for i in range(18)]


# start preprocessing
count = 0
expected_count = len(phi_list)*len(Ca_list)*len(angle_list)
progress = 0

f_main = open("main.log", "a+")
f_pre = open("preprocess.log", "a+")

for phi in phi_list:
    for Ca in Ca_list:
        for angle in angle_list:
            fn = "phi_{}_Re_0.1_Ca_{}_aggregation_1KT_ncycle_{}_np_2_angle_{}".format(phi, Ca, nstrain, angle)

            subprocess.call(["python3", "Preprocessing.py", fn, "system=0", "verbose=0"], stdout = f_pre, stderr = f_pre)
            count += 1

            if count/expected_count >= progress:
                f_main.write("Progress: {} % ({}/{}) at {}\n".format(round(progress*100), count, expected_count, datetime.datetime.now()))
                f_main.flush()
                progress += 0.05

f_main.write("Finish!\n")
f_main.write("Total time elapsed = {}\n".format(str(datetime.timedelta(seconds=time.time()-start_time))))
f_main.flush()
f_main.close()
f_pre.close()
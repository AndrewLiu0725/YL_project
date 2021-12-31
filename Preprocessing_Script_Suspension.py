# ===============================================================================
# Copyright 2021 An-Jun Liu
# Last Modified Date: 12/29/2021
# ===============================================================================
import os
import subprocess
import time
import datetime

'''
This program is a script for preprocessing the data for suspension system
'''

start_time = time.time()

path = "/raid6/ctliao/Data/HI_ordering/"

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


# collect the job name
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


phis = []

number_of_parameter_set = 0
for phi in parameter_set.keys():
    if (len(parameter_set[phi].keys()) > 5):
        for Ca in parameter_set[phi].keys():
            number_of_parameter_set += len(parameter_set[phi][Ca])
        phis.append(phi)

# Preprocessing
count = 0
expected_count = number_of_parameter_set
progress = 0

f_main = open("main.log", "a+")
f_pre = open("preprocess.log", "a+")

for fn in os.listdir(path):
    if (fn[0] == "h") and (os.path.isdir(path+fn)):
        result = parser(fn)
        if result[0] in phis:
            subprocess.call(["python3", "Preprocessing.py", fn, "system=1", "verbose=0", "check=0"], stdout = f_pre, stderr = f_pre)
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
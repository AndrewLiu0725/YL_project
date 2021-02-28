# ===============================================================================
# Copyright 2021 An-Jun Liu
# Last Modified Date: 02/28/2021
# ===============================================================================
import subprocess
import sys
import time
import datetime
from CheckJob import checkJobs

"""
This code is the script for sending jobs to the server.
It will send the progress and check whether the submitted jobs are complete every 5% jobs are submitted.
This script also ensures that it will not overwhelm the server by limiting the idle jobs under 36.
"""

# setup
nstrain = 4000
vol = 746.3163
test = 0 # 1 for creating init case for each set of parameters
Ca_list = [round(i*0.01+0.01, 5) for i in range(20)]
phi_list = [2*vol/(24*x**2)*100 for x in range(30, 41)] # y is fixed to 24
angle_list = [90 - 10*i for i in range(18)]
queue_size = 36 # number of idle jobs
count = 0
expected_count = len(phi_list)*len(Ca_list)*len(angle_list)

progress = 0

if not test:
	f = open("main.log", "a+")
	for phi in phi_list:
		for Ca in Ca_list:
			for angle in angle_list:
				stat = subprocess.check_output(["condor_q", "ajliu"]).decode('utf-8')
				idle_count = int(stat.split('\n')[-2].split()[6])

				# busy waiting for entering the queue
				while idle_count >= queue_size:
					time.sleep(1*60) # wait for 1 min
					stat = subprocess.check_output(["condor_q", "ajliu"]).decode('utf-8')
					idle_count = int(stat.split('\n')[-2].split()[6])

				# send job	
				subprocess.call(['python3', 'Parameter.py', 'phi_{}_Re_0.1_Ca_{}_aggregation_1KT_ncycle_{}_np_2_angle_{}'.format(phi, Ca, nstrain, angle)])
				count += 1

				# print progress and check jobs every 5% jobs are submitted
				if count/expected_count >= progress:
					#print("Progress: {} % ({}/{})".format(progress*100, count, expected_count), flush = True)
					f.write("Progress: {} % ({}/{})\n".format(int(progress*100), count, expected_count))
					f.flush()
					checkJobs()
					progress += 0.05
	
	time.sleep(5*60)
	checkJobs()
	f.write("Finish!\n")
	f.flush()
	f.close()
		

else:
	for x in range(30, 41):
		if not(x in [36, 40]): continue
		phi = 2*vol/(24*x**2)*100
		subprocess.call(['python3', 'Parameter.py', 'test_phi_{}_Re_0.1_Ca_{}_aggregation_1KT_ncycle_{}_np_2'.format(phi, 0.01, 1)])
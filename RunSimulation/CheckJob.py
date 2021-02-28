# ===============================================================================
# Copyright 2021 An-Jun Liu
# Last Modified Date: 02/28/2021
# ===============================================================================
import os
import subprocess
import datetime

"""
This function is to check whether the submitted jobs are run successfully.
If incomplete or failed, resubmit it again.
"""

def checkJobs():
	f = open("main.log", "a+")
	f.write('--------------------------------------------------------------------------\n')
	# get the id of current running jobs to prevent delete their data
	running_jobs = subprocess.check_output(["condor_q", "ajliu"]).decode('utf-8')
	running_jobs_id = []
	for line_index, line in enumerate(running_jobs.split('\n')):
		if line_index > 3:
			if len(line.split()) == 0:
				break
			running_jobs_id.append(int(float(line.split()[0])))

	# recognize the incomplete or failed jobs and resubmit them
	root_folder = "/userdata4/ajliu/RBC_doublet/"
	resubmit_count = 0

	for fn in os.listdir(root_folder+"Data"):
		if (fn.split('_')[0] != 'phi'):
			continue

		log_id = [] # list of log id for this case

		for subfn in os.listdir(root_folder+'Data/'+fn):
			if subfn.split('.')[0] == 'log':
				log_id.append(int(subfn.split('.')[1]))

		if max(log_id) in running_jobs_id: # skip if it is running right now
			continue

		job_completed = 0

		# chech whether this job run successfully
		ncycle = int(fn.split('_')[9])
		if os.path.isfile(root_folder+"Data/"+fn+"/data/bond0_t{}.vtk".format(int(3669*2*(int(ncycle/2)-1)))):
			job_completed = 1

		if not job_completed:
			#print("Resubmit job {} at {}".format(fn, datetime.datetime.now()), flush = True)
			f.write("Resubmit job {} at {}\n".format(fn, datetime.datetime.now()))
			f.flush()
			subprocess.call(['rm', '-r', root_folder+'Data/'+fn])
			com = fn.split('_')
			com[1] = str(float(com[1]))	
			subprocess.call(["python3", root_folder+"Script/Parameter.py", '_'.join(com)])
			resubmit_count += 1
	#print('\n--------------------------------------------------------------------------', flush = True)
	#print("Resubmit {} {} at {}".format(resubmit_count, "jobs" if resubmit_count > 1 else "job",datetime.datetime.now()), flush = True)
	f.write("Resubmit {} {} at {}\n".format(resubmit_count, "jobs" if resubmit_count > 1 else "job",datetime.datetime.now()))
	f.write('--------------------------------------------------------------------------\n')
	f.flush()
	f.close()


if __name__ == "__main__":
    checkJobs()
import sys
import subprocess
import os

# Create .dat file
subprocess.call(['python3', 'parameter.py', sys.argv[1]])

components = sys.argv[1].split('_')

for fn in os.listdir("./Data/"):
    if fn.split('_')[0:4] == sys.argv[1].split('_')[0:4] and int(sys.argv[1].split('_')[-1]) > int(fn.split('_')[-1]):
        cycle_number = int(sys.argv[1].split('_')[-1]) - int(fn.split('_')[-1])
        folder_to_copy = fn

# Scripts
if len(components) == 6:
    try:
        subprocess.call(['mv', 'Data/'+folder_to_copy, 'Data/'+sys.argv[1]])
        print('Complete copying the config of particles from {}'.format(folder_to_copy))
    except:
        folder_particle = '_'.join(sys.argv[1].split('_')[0:2]) + '_ncycle_2'
        folder_new = sys.argv[1]
        subprocess.call(['cp','-r','Data/'+folder_particle, 'Data/'+folder_new])
        print('Complete copying the config of particles from {}'.format(folder_particle))
print('Start running the simulation "{}"'.format(sys.argv[1]))
subprocess.call(['sh', 'job.sh', sys.argv[1]])                                                                                                                               

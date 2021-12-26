import subprocess

filename = "phi_6.9_Re_0.1_Ca_0.05_aggregation_1KT_ncycle_4000_np_2_angle_10"
filetype = ["parameter.txt", "COMs.npy", "COMs_NB.npy", "Ypos_t.npy"]

for i in range(4):
    subprocess.call(['cp', '/Users/andrewliu/remote_disk/Data_Transfer/{}_{}'.format(filename, filetype[i]), './data/'])
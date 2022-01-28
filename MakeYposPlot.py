# ===============================================================================
# Copyright 2021 An-Jun Liu
# Last Modified Date: 10/17/2021
# ===============================================================================
import numpy as np 
import matplotlib.pyplot as plt

SHOW = 0 # 0 for save figures, 1 for show plots
plt.rcParams['font.family'] = 'DeJavu Serif'
plt.rcParams['font.serif'] = ['Times New Roman']

data = np.load("./Data/phi_4.7_Re_0.1_Ca_0.01_aggregation_1KT_ncycle_4000_np_2_angle_0_Ypos_t.npy")
#data = np.load("../../remote_disk/Data_Transfer/phi_6.9_Re_0.1_Ca_0.2_aggregation_1KT_ncycle_4000_np_2_angle_60_Ypos_t.npy")

fig, ax = plt.subplots(figsize = (9, 6))
st, et = 50, 100
indices = [7,8,11]
for i in range(3):
    ax.plot((np.arange(2000)*2)[st:et], data[:, indices[i]][st:et]*0.5, label = "node {}".format(i+1))
ax.set_xlabel("{}".format(r'$\dot \gamma t$'), fontsize = 20)
ax.set_ylabel("y ({})".format(r'$\mu m$'), fontsize = 20)
ax.legend(prop={'size': 15}, loc = 1, frameon = False)
#ax.set_title("Selected Nodes' y(t)", fontsize = 20)
ax.xaxis.set_tick_params(labelsize = 15)
ax.yaxis.set_tick_params(labelsize = 15)
fig.tight_layout()

if SHOW:
    plt.show()
else:
    plt.savefig("./Pictures/Manuscript/Manuscript_YPos_1.png", dpi = 300)


fig, ax = plt.subplots(figsize = (9, 6))
offset = int(0.0135*4000 + 1) # igonre the physically impossible high frequency regime
fs = 0.5 
N = 2000
for i in range(3):
    yf = np.fft.rfft(data[:, indices[i]] - np.mean(data[:, indices[i]]))
    xf = np.fft.rfftfreq(N, d = 1.0/fs)
    ax.plot(xf[offset:], (2.0/N * np.abs(yf))[offset:], label = "node {}".format(i+1))
ax.set_xlabel("frequency ({})".format(r'$\dfrac{1}{\dot \gamma t}$'), fontsize = 20)
ax.set_ylabel("Amplitude", fontsize = 20)
ax.xaxis.set_tick_params(labelsize = 15)
ax.yaxis.set_tick_params(labelsize = 15)
ax.legend(prop={'size': 15}, loc = 1, frameon = False)
#ax2.set_title("Corresponding Frequency Spectrum", fontsize = 20)
fig.tight_layout()
if SHOW:
    plt.show()
else:
    plt.savefig("./Pictures/Manuscript/Manuscript_YPos_2.png", dpi = 300)
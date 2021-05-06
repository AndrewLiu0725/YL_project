# ===============================================================================
# Copyright 2021 An-Jun Liu
# Last Modified Date: 04/16/2021
# ===============================================================================
import numpy as np 
import matplotlib.pyplot as plt

data = np.load("./Data/phi_6.0_Re_0.1_Ca_0.01_aggregation_1KT_ncycle_4000_np_2_angle_30_Ypos_t.npy")

fig, (ax1, ax2) = plt.subplots(1, 2, figsize = (16, 7))
st, et = 50, 100
indices = [7,8,11]
for i in range(3):
    ax1.plot((np.arange(2000)*2)[st:et], data[:, indices[i]][st:et]*0.5, label = "node {}".format(i+1))
ax1.set_xlabel("time ({})".format(r'$\dot \gamma t$'), fontsize = 20)
ax1.set_ylabel("y ({})".format(r'$\mu m$'), fontsize = 20)
ax1.legend(prop={'size': 12})
ax1.set_title("Selected Nodes' y(t)", fontsize = 20)
ax1.xaxis.set_tick_params(labelsize = 12)
ax1.yaxis.set_tick_params(labelsize = 12)


offset = 10 # igonre the physically impossible high frequency regime
fs = 0.5 
N = 2000
for i in range(3):
    yf = np.fft.rfft(data[:, indices[i]])
    xf = np.fft.rfftfreq(N, d = 1.0/fs)
    ax2.plot(xf[offset:], (2.0/N * np.abs(yf))[offset:], label = "node {}".format(i+1))
ax2.set_xlabel("frequency ({})".format(r'$\dfrac{1}{\dot \gamma t}$'), fontsize = 20)
ax2.set_ylabel("Amplitude", fontsize = 20)
ax2.xaxis.set_tick_params(labelsize = 12)
ax2.yaxis.set_tick_params(labelsize = 12)
ax2.legend(prop={'size': 12})
plt.title("Corresponding Frequency Spectrum", fontsize = 20)
fig.tight_layout()
#plt.show()
plt.savefig("./Pictures/Manuscript_YPos.png", dpi = 300)
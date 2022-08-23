# ===============================================================================
# Copyright 2021 An-Jun Liu
# Last Modified Date: 08/22/2022
# ===============================================================================
import numpy as np 
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

SHOW = 0 # 0 for save figures, 1 for show plots
plt.rcParams['font.family'] = 'DeJavu Serif'
plt.rcParams['font.serif'] = ['Times New Roman']

data = np.load("./Data/phi_4.7_Re_0.1_Ca_0.01_aggregation_1KT_ncycle_4000_np_2_angle_0_Ypos_t.npy")
#data = np.load("../../remote_disk/Data_Transfer/phi_6.9_Re_0.1_Ca_0.2_aggregation_1KT_ncycle_4000_np_2_angle_60_Ypos_t.npy")

rows, columns = 1, 2
fig, axs = plt.subplots(rows, columns, figsize = (9*columns, 6*rows))

# left plot
ax = axs[0]
st, et = 50, 100
indices = [7,8,11]
for i in range(3):
    ax.plot((np.arange(2000)*2)[st:et], data[:, indices[i]][st:et]*0.5, label = "node {}".format(i+1))
ax.set_xlabel("{}".format(r'$\gamma _{rot}$'), fontsize = 20)
ax.set_ylabel("y ({})".format(r'$\mu m$'), fontsize = 20)
#ax.set_title("Selected Nodes' y(t)", fontsize = 20)
ax.tick_params(which='both', labelsize=15, width=2, length=8, direction='in')
ax.xaxis.set_major_locator(MaxNLocator(5))
ax.xaxis.set_minor_locator(MaxNLocator(10))
ax.yaxis.set_major_locator(MaxNLocator(5)) 
#ax.yaxis.set_minor_locator(MaxNLocator(10))
ax.legend(prop={'size': 15}, loc = 1, frameon = False)


# right plot
ax = axs[1]
offset = int(0.0135*4000 + 1) # igonre the physically impossible high frequency regime
fs = 0.5 
N = 2000
for i in range(3):
    yf = np.fft.rfft(data[:, indices[i]] - np.mean(data[:, indices[i]]))
    xf = np.fft.rfftfreq(N, d = 1.0/fs)
    ax.plot(xf[offset:], (2.0/N * np.abs(yf))[offset:], label = "node {}".format(i+1))
ax.set_xlabel("frequency ({})".format(r'$\dfrac{1}{\gamma _{rot}}$'), fontsize = 20)
ax.set_ylabel("Amplitude", fontsize = 20)
ax.tick_params(which='both', labelsize=15, width=2, length=8, direction='in')
ax.xaxis.set_major_locator(MaxNLocator(5))
ax.xaxis.set_minor_locator(MaxNLocator(10))
ax.yaxis.set_major_locator(MaxNLocator(5)) 
#ax.yaxis.set_minor_locator(MaxNLocator(10))
ax.legend(prop={'size': 15}, loc = 1, frameon = False)
#ax2.set_title("Corresponding Frequency Spectrum", fontsize = 20)

fig.tight_layout()
if SHOW:
    plt.show()
else:
    fig.savefig("./Pictures/Manuscript/fig_S3.png", dpi = 300)
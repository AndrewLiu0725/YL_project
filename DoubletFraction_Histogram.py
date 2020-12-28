from Doublet_Functions import calcDoubletFraction
import matplotlib.pyplot as plt
from scipy import stats
from sklearn.metrics import mean_squared_error
import numpy as np 

"""
This code is to plot doublet fraction's histogram and fit it with skew normal distribution in suspension system
"""


phi, Ca = 2.3994, 0.01
result = calcDoubletFraction(phi, Ca, 1, 1, 0)
data = (result[0][0]*result[4]/2).astype(int) # convert to number of doublets
fig = plt.figure(figsize=(8, 6))
ax1 = fig.add_subplot(111)    
ax2 = ax1.twiny()
max_num_of_doublet = max(data)
min_num_of_doublet = min(data)
hist_data = ax1.hist(data, bins = (max_num_of_doublet-min_num_of_doublet+1), density = True, label = "Simulation data") # plot pmf
if (len(hist_data[0].nonzero()) == 1):
    pass
else:
    # fit with skew normal
    a, loc, scale = stats.skewnorm.fit(data)
    x = np.linspace(min_num_of_doublet, max_num_of_doublet, 100)
    p = stats.skewnorm.pdf(x, a, loc, scale)
    ax1.plot(x, p, 'k', linewidth=2, label = "Fitted line")
    # calculate rms error
    x_predicted = np.array(list(range(min_num_of_doublet, max_num_of_doublet + 1)))
    df_estimated = stats.skewnorm.pdf(x_predicted, a, loc, scale)
    rms = np.sqrt(mean_squared_error(hist_data[0], df_estimated))
    ax1.text(max_num_of_doublet*0.85, max(hist_data[0])*0.6, "RMSD = {:.3}\na = {:.3}\nloc = {:.3}\nscale = {:.3}".format(rms, a, loc, scale), weight = "bold")

print(sum(hist_data[0]))

# set second x axes
def tick_function(x):
    return 2*x/result[4]

ax2.set_xbound(tick_function(ax1.get_xbound()[0]), tick_function(ax1.get_xbound()[1]))
ax1.set_xlabel("# of doublets", fontsize = 15) 
ax1.set_ylabel("Frequency", fontsize = 15)
ax1.set_title("Doublet fraction's histogram\nphi = {}, Ca = {} (whole time series)".format(phi, Ca), fontsize = 20)
ax1.legend()
ax2.set_xlabel('Doublet fraction', fontsize = 15)
fig.tight_layout()
plt.show()
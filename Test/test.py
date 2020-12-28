from scipy import stats
import matplotlib.pyplot as plt
import numpy as np

# choose some parameters
a, loc, scale = 2, 0, 1
# draw a sample
data = stats.skewnorm(a, loc, scale).rvs(1000)
# estimate parameters from sample
ae, loce, scalee = stats.skewnorm.fit(data)
# Plot the PDF.
plt.figure()
plt.hist(data, bins=100, density=True, alpha=0.6, color='g')
xmin, xmax = plt.xlim()
x = np.linspace(xmin, xmax, 100)
p = stats.skewnorm.pdf(x,ae, loce, scalee)#.rvs(100)
plt.plot(x, p, 'k', linewidth=2)
plt.show()
import numpy as np 
import matplotlib.pyplot as plt

"""
This code is to plot eq 37 in the following paper: 
https://smartech.gatech.edu/bitstream/handle/1853/52903/Faroughi_Huber_Accepted_manuscript.pdf?sequence=1&isAllowed=y
"""

def dilute_formula(Ca, phi):
    l = 1
    k = (((2*l+3)*(19*l+16))/(40*(l + 1)))**2
    return (1/(1+k*(Ca**2)))*(phi/(1-phi))*(((1+2.5*l)/(1+l))+((140*(l**3+l**2-l-1))/(28*(2*l+3)*((l+1)**2)))*k*(Ca**2))

nop = 20
Cas = np.linspace(0.01, 0.2, nop)
iv = np.zeros(nop)
for i in range(nop):
    iv[i] = dilute_formula(Cas[i], 0.05)
plt.plot(Cas, iv, "o")
plt.xlabel("Ca")
plt.ylabel(r'$\left[ \eta \right]$')
plt.show()
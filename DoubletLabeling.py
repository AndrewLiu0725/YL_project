# ===============================================================================
# Copyright 2021 An-Jun Liu
# Last Modified Date: 01/10/2022
# ===============================================================================
from RBC_Utilities import calcDoubletFraction
import numpy as np 

"""
This code creates the vtk file with doublets only (need the help of RBC_Utilities.py's doublet_or_not ndarray)
Users can superpose the output vtk files (doublet0_{}.vtk) on the original vtk files to see the doublet counting results.
"""

# phi_range [2.3994, 2.9993, 3.4492, 3.8991, 4.9488, 5.9986]
phi, Ca, ensemble_id = 3.4492, 0.14, 6
filepath = 'data/h24_phi{}_Re0.1_Ca{}_WCA1_zero0.8-{}/'.format(phi, Ca, ensemble_id)
df_result = calcDoubletFraction(phi, Ca, 1, 1, 0, ensemble_id, 1, outputDataType=1)

for t in range(1825, 1975):
    doublet_indices = df_result[2][np.where(df_result[1][:, t])].flatten()

    with open(filepath+'particle0_{}.vtk'.format(t*4000), 'r') as f:
        data = f.readlines()

    numParticles = int(int(data[5].split()[1])/642)
    numDoublets = len(doublet_indices)
    
    with open(filepath+'doublet0_{}.vtk'.format(t*4000), 'w') as f:
        f.writelines(data[:5])

        f.write('POINTS {} float\n'.format(642*numDoublets))

        for ptcl in doublet_indices:
            f.writelines(data[6+ptcl*642: 6+(ptcl+1)*642])
    
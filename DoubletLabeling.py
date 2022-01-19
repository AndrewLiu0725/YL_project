# ===============================================================================
# Copyright 2021 An-Jun Liu
# Last Modified Date: 01/19/2022
# ===============================================================================
from RBC_Utilities import calcDoubletFraction
import numpy as np 

"""
This code creates the vtk file with doublets only (need the help of RBC_Utilities.py's doublet_or_not ndarray)
Users can superpose the output vtk files (doublet0_{}.vtk) on the original vtk files to see the doublet counting results.
"""

def labelDoublet(phi, Ca, dependent, system, st, et):
    '''
    system = 0 means two-cell system and 1 mens suspension 
    st and et are in vtk file unit
    '''
    if system:
        filepath = 'data/h24_phi{}_Re0.1_Ca{}_WCA1_zero0.8-{}/'.format(phi, Ca, dependent)
        df_result = calcDoubletFraction(phi, Ca, 1, 1, 0, dependent, 1, outputDataType=1)
    else:
        filepath = 'data/phi_{}_Re_0.1_Ca_{}_aggregation_1KT_ncycle_4000_np_2_angle_{}/data/'.format(phi, Ca, dependent)
        df_result = calcDoubletFraction(phi, Ca, 1, 1, 0, dependent, 0, outputDataType=1)

    interval = 4000 if system else 7338
    
    for t in range(st, et):
        doublet_indices = df_result[2][np.where(df_result[1][:, t if system else 2*t])].flatten()

        vtk_t = t*interval if t else '0000'
        with open(filepath+'{}{}.vtk'.format('particle0_' if system else 'bond0_t', vtk_t), 'r') as f:
            data = f.readlines()

        #numParticles = int(int(data[5].split()[1])/642)
        numDoublets = len(doublet_indices)
        
        with open(filepath+'doublet0_{}.vtk'.format(vtk_t), 'w') as f:
            f.writelines(data[:5])

            f.write('POINTS {} float\n'.format(642*numDoublets))

            for ptcl in doublet_indices:
                f.writelines(data[6+ptcl*642: 6+(ptcl+1)*642])

if __name__ == "__main__":
    labelDoublet(3.8, 0.04, -70, 0, 0, 1999)
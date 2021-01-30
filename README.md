# RBC Project Github Repo
## Introduction
This project is An-Jun Liu's work in the Polymer Physics and Complex Fluids Group led by Professor Yeng-Long Chen at Institute of Physics, Academia Sinica.
The two main contributions of this project is:
1. Efficiently calculate the time series of the number of doublets in a fluid system. The utility function also outputs the time series of states (doulet state, two-singlet state, and sliding state) of a fluid system. Before I wrote this utility function, our group members could only know the doublet fraction by visualizing the fluid system and watch the animationm, which is very time-consuming.
2. Doing data analysis based on the utility function mentioned above. The analyses include the relationship between intrinsic/relative viscosity and volume fraction/capilary number, the relationship between doublet fraction and intrinsic viscosity, and the relationship between elastic stress tensor and interparticle distance (i.e. a charatersitic of doublet formation).

## Prerequisite
To use the utility function in this project, you have to do the following:
1. Install python3
2. Install numpy
3. Install scipy
4. Run this command "gcc -fPIC -shared -o RBC_Utilities_CExtension.so RBC_Utilities_CExtension.c" to compile and create .so file in the folder where you clone this repo

## Code descriptions
### RBC_Utilities.py
This utility python file provides functions to calculate doublet fraction (accelerated by preprocessing and writing a C extension file), interparticle/elastic stress, intrinsic viscosity, and relative viscosity for two-cell and suspension system.

To use calcDoubletFraction(), you have to preprocess the data by using Preprocess.py.

### Preprocess.py
This file is to preprocess the data for calculating the doublet function.
The preprocessed data includes
1. time series of bounded center of mass of each particle, i.e. _COMs.npy
2. time series of unbounded center of mass of each particle, i.e. _COMs_NB.npy
3. time series of y coordinate of specific nodes, i.e. _Ypos_t.npy
4. parameters including timesteps (COMs time unit), interval (Ypos time unit), particle numbers, points per particle (used in Ypos_t), and dimensions, i.e. _parameter.txt

### Data analysis
#### DoubletFraction_Deviation.py
Calculate the deviation (in time, after 500 strains) of each ensembled average doublet fraction time series
#### DoubletFraction_Histogram.py
Plot doublet fraction's histogram and fit it with skew normal distribution in suspension system
#### DoubletFraction_vs_Time_Suspension_EnsembleAveraged.py
Plot ensemble averaged doublet fraction time series for suspension system
#### ElasticStress_InterParticleDistance_DoubletFraction_vs_Time.py
Plot elastic stress tensor, interpartilce distance, and doublet fraction/state time series
#### ElasticStress_vs_InterParticleDistance.py
Plot elastic stress tensor vs interparticle distance (avg or max).
Also provides linear regression line
#### PhaseDiagram_Suspension_EnsembleAveraged.py
Make two plots: doublet fraction vs Ca plot and its phasediagram (in colormap)
#### PhaseDiagram_TwoCell.py
Plot phase diagram for two-cell system.
#### Slope_vs_Ca_Suspension_EnsembleAveraged.py
Plot slope (defined by doublet fraction vs intrinsic viscosity) vs Ca for suspension system.
#### Slope_vs_Ca_TwoCell.py
Plot slope (defined by doublet fraction vs intrinsic viscosity) vs Ca for two-cell system.
#### Viscosity_vs_Ca_SeperatedState.py
Plot intrinsic and relative viscosity vs Ca for both doublet and two-singlet state in two-cell system with all volume fraction.
Note that the states are seperated!
#### Viscosity_vs_Ca_UnseperatedState.py
Plot intrinsic and relative viscosity vs Ca for two-cell system with all volume fraction.
Note that states are not seperated!
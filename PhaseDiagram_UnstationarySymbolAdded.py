# ===============================================================================
# Copyright 2021 An-Jun Liu
# Last Modified Date: 08/16/2022
# ===============================================================================
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from matplotlib import cm
from matplotlib.ticker import MultipleLocator
import pickle

"""
This code is to add unstationary symbol on the phase diagram
"""

DEBUG = 0

SHOW = 0
print("Setting:\nSHOW:{}".format(SHOW))

# utility functions
def Z(x, y, df):
    return df[x, y]

def findBorder(target, value):
    left = 0
    right = len(target)

    reverse_order = 1 if target[0] > target[1] else 0 # 0 for ascending, 1 for descending

    if reverse_order:
        while True:
            mid = int((left+right)/2)
            if target[mid] >= value and target[mid+1] <= value:
                break
            elif target[mid] <= value:
                right = mid
            else:
                left = mid+1

    else:
        while True:
            mid = int((left+right)/2)
            if target[mid] <= value and target[mid+1] >= value:
                break
            elif target[mid] >= value:
                right = mid
            else:
                left = mid+1
    return [target[mid], target[mid+1]]

def makePlot(variable_type, system, division_type, unstationary_symbol_style, alpha):
    """
    Input:

    unstationary_symbol_style: 1 means same-size marker, 0 means marker which occupy the whole rectangle
    """

    # setup
    threshold = 0.1 # threshold for the ratio of the set of parameters to be marked as unstationry
    print("Current task:\nvariable_type:{}, system:{}, division_type:{}, unstationary_symbol_style:{}, alpha:{}"
    .format(variable_type, system, division_type, "equal" if unstationary_symbol_style else "filled", alpha))



    # read the data for phase diagram (grid phi, grid Ca, df)
    with open("Data/PhaseDiagram_{}_new.pickle".format(system), 'rb') as handle:
        data_PD = pickle.load(handle)

    # get the positions of the unstationary cases
    with open("Data/UnstationaryRatio{}{}_{}_alpha_{}_p_20.pickle".format("_IV_" if variable_type == "IV" else "_", system, division_type, alpha), 'rb') as handle:
        unstationary_ratio = pickle.load(handle) # format: [phi, Ca, unstationary ratio]
        unstationary_case = []
        for point in unstationary_ratio:
            if point[2] > threshold: # mark as unstationary
                unstationary_case.append([point[0], point[1]])

    ### Plot
    # plot the phase diagram
    plt.rcParams['font.family'] = 'DeJavu Serif'
    plt.rcParams['font.serif'] = ['Times New Roman']
    fig, ax = plt.subplots(figsize = (16,12))
    [grid_phi, grid_Ca, df, index_X, index_Y] = data_PD
            
    X, Y = np.meshgrid(grid_phi, grid_Ca)
    im = ax.pcolormesh(Y, X, Z(index_X, index_Y, df), cmap=cm.coolwarm, vmin=0, vmax=0.6 if system == "Suspension" else 0.6)
    cb = fig.colorbar(im, ax = ax)
    cb.ax.tick_params(labelsize = 25)
    #ax.set_title("Phase Diagram of {} vs Ca ({} system)".format(r'$\phi$', system), fontsize = 30)

    # add unstationary symbol
    if unstationary_symbol_style:
        for point in unstationary_case:
            ax.plot(point[1], point[0], marker = "x", color = "k", markersize = 20)
    else:
        for point in unstationary_case: # [phi, Ca]
            Ca_border = findBorder(grid_Ca, point[1])
            phi_border = findBorder(grid_phi, point[0])
            ax.plot([Ca_border[0], Ca_border[1]], [phi_border[0], phi_border[1]], "k")
            ax.plot([Ca_border[0], Ca_border[1]], [phi_border[1], phi_border[0]], "k")

    void_marker = mlines.Line2D([], [], color = "k", marker = "x", linestyle = 'None', markersize = 20, label="unstationary")

    ax.set_xlabel('Ca', fontsize = 35)
    ax.set_ylabel(r'$\phi$'+"(%)", fontsize = 40)
    ax.tick_params(which='both', labelsize=30, width=4, length=12, direction='in')
    ax.xaxis.set_major_locator(MultipleLocator(0.05))
    ax.xaxis.set_minor_locator(MultipleLocator(0.025))
    ax.yaxis.set_major_locator(MultipleLocator(1))
    ax.yaxis.set_minor_locator(MultipleLocator(0.5))
    ax.legend(handles = [void_marker], bbox_to_anchor=(0.9, 1.07), loc='upper left', prop={'size': 25}, frameon = False)

    fig.tight_layout()
    
    if SHOW:
        plt.show()
    else:
        plt.savefig('./Pictures/PhaseDiagram_UnstationarySymbolAdded/{}System_PhaseDiagram_UnstationarySymbolAdded_{}_{}_alpha_{}_SymbolStyle_{}.png'.format(
        system, "IntrinsicViscosity" if variable_type == "IV" else "DoubletFraction",division_type, alpha, unstationary_symbol_style), dpi = 200)



if __name__ == "__main__":
    system_list = ["Suspension", "TwoCell"]
    #system_list = ["Suspension"]
    #division_type_list = ["EnsembleAveraged", "Indivisual"]
    division_type_list = ["EnsembleAveraged"]
    for system in system_list:
        for division_type in division_type_list:
            for alpha in [0.1]:
                for unstationary_symbol_style in [0]:
                    makePlot("IV", system, division_type, unstationary_symbol_style, alpha)
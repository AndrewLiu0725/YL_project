# ===============================================================================
# Copyright 2021 An-Jun Liu
# Last Modified Date: 01/30/2021
# ===============================================================================
import numpy as np 
import matplotlib.pyplot as plt
from RBC_Utilities import calcDoubletFraction, getInstrinsicViscosity, getStress

"""
This code is to plot elastic stress tensor, interpartilce distance, and doublet fraction/state time series.
"""

def make_patch_spines_invisible(ax):
    ax.set_frame_on(True)
    ax.patch.set_visible(False)
    for sp in ax.spines.values():
        sp.set_visible(False)

def plot(plot_or_not):
    fig, host = plt.subplots(figsize = (12, 9))
    fig.subplots_adjust(right=0.75)

    par1 = host.twinx()
    par2 = host.twinx()
    par2.spines["right"].set_position(("axes", 1.1))
    make_patch_spines_invisible(par2)
    par2.spines["right"].set_visible(True)

    result = calcDoubletFraction(phi, Ca, 1.0, r, 0, angle, 0)
    dfg = result[0][0][:ncycle]
    simul_time = list(range(len(dfg)))
    #p1, = host.plot(simul_time[st : et], dfg[st : et], "b-", label = "doublet fraction")
    p1, = host.plot(simul_time[st : et], result[2][0][0, st : et], "b-", label = "state series")
    p2, = par1.plot(simul_time[st : et], result[3][0,st : et], "r-", label = "interparticle distance")
    #p1, = host.plot(simul_time[st : et], result[6][0,st : et, 0], "r-", label = "particle x position", ls = '--')
    p3, = par2.plot(simul_time[st : et], getStress(phi, Ca, 0, 0, ncycle, angle, 0)[:,3][st : et], "g-", label = r'$\sigma _{elas, yx}$')

    host.set_title('phi = {}, Ca = {}, angle = {}\ncriteria r = {}Dm'.format(phi, Ca, angle, r), fontsize = 30)
    host.set_xlabel("timesteps(strain)", fontsize = 20)
    host.set_ylabel("Doublet fraction", fontsize = 20)
    #host.set_ylim([-0.1, 1.1])
    par1.set_ylabel("interparticle distance", fontsize = 20)
    #par1.set_ylabel(r'$\left[ \eta \right]$', fontsize = 20)
    par2.set_ylabel(r'$\sigma _{elas, yx}$', fontsize = 20)

    host.yaxis.label.set_color(p1.get_color())
    par1.yaxis.label.set_color(p2.get_color())
    par2.yaxis.label.set_color(p3.get_color())

    tkw = dict(size=6, width=1.5)
    host.tick_params(axis='y', colors=p1.get_color(), **tkw)
    par1.tick_params(axis='y', colors=p2.get_color(), **tkw)
    par2.tick_params(axis='y', colors=p3.get_color(), **tkw)
    host.tick_params(axis='x', **tkw)

    lines = [p1, p2, p3]
    host.legend(lines, [l.get_label() for l in lines], prop={'size': 12}, bbox_to_anchor=(1.2,1.2))
    fig.tight_layout()
    if plot_or_not == 1:
        plt.savefig("./Pictures/TwoCellSystem_TransionalPart_stress_yx_phi_{}_Ca_{}_angle_{}_r_{}_start_{}_end_{}.png".format(phi, Ca, angle, r, st, et), dpi = 300)
    plt.show()

r = 1
ncycle = 2000
phi, Ca, angle, st, et = 7.9, 0.04, -40, 0, 1900
plot(0)

'''
3.8, 0.03, -40, 0, 200
4.7, 0.04, -70, 100, 300
6.9, 0.04, 80, 1050, 1400
6.9, 0.05, 10, 150, 350
7.9, 0.03, 20, 550, 850
7.9, 0.04, 70, 250, 450
7.9, 0.05, -40, 550, 750
6.9, 0.06, -40, 700, 900
6.9, 0.06, -10, 50, 350
6.0, 0.06, 80, 600, 850
6.0, 0.06, -80, 150, 500
6.0, 0.05, 10, 50, 300
4.7, 0.06, 40, 1250, 1500
3.8, 0.05, 30, 0, 250
'''
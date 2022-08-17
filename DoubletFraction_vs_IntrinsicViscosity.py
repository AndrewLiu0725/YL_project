# ===============================================================================
# Copyright 2021 An-Jun Liu
# Last Modified Date: 08/16/2022
# ===============================================================================
#from tkinter import font
import matplotlib.pyplot as plt
import numpy as np 
from RBC_Utilities import calcDoubletFraction, getIntrinsicViscosity, getSuspensionParameterSets
from matplotlib.ticker import MaxNLocator

def makeDFvsIVplot(ps, plotType, plotTitle, plotName):
    '''
    plotType: 0 for making scatter plots, 1 for averaging intrinsic viscosity w.r.t. to individual doublet fraction before making the plot
    '''
    
    #plt.rcParams['font.family'] = 'DeJavu Serif'
    #plt.rcParams['font.serif'] = ['Times New Roman']
    
    fig, ax = plt.subplots(figsize = (8, 6))
    [_, parameter_set] = getSuspensionParameterSets()
    st = 5

    for [phi, Ca] in ps:
        if len(plotTitle) > 0:
            parameter = '{}={:.2}%'.format(r'$\phi$', phi) if plotTitle[0] == 'C' else 'Ca={}'.format(Ca)
        else:
            parameter = ''

        if plotType == 0:
            df_total, iv_total = [], [] # for making linear regression
        else:
            data = {}           

        for ensemble_id in parameter_set[phi][Ca]:
            try:
                [df, et, _] = calcDoubletFraction(phi, Ca, 1, 1, 0, ensemble_id, 1)
                iv = getIntrinsicViscosity(phi, Ca, None, ensemble_id, 1)

                if plotType == 0:
                    # in this case, the linear regression is made to the whole data points
                    ax.scatter(df[0][st:et], iv[st:et], color='b', facecolor='none')
                    df_total += list(df[0][st:et])
                    iv_total += list(iv[st:et])
                else:
                    # sort out iv[t] according to df[t]
                    for t in range(st, et):
                        if df[0][t] in data.keys():
                            data[df[0][t]].append(iv[t])
                        else:
                            data[df[0][t]] = [iv[t]]
            except:
                pass

        if plotType == 1:
            df_total = list(data.keys())
            df_total.sort()
            iv_total = [np.mean(data[df]) for df in df_total]
            iv_std = [np.std(data[df]) for df in df_total]
            #ax.errorbar(df_total, iv_total, iv_std, color='b', capsize = 2, marker = 's')
            ax.plot(df_total, iv_total, marker='s', markerfacecolor='None', linestyle='None')

        model = np.polyfit(df_total, iv_total, 1)
        predict = np.poly1d(model)
        x_fit = np.linspace(min(df_total), max(df_total), 100)
        y_fit = predict(x_fit)
        ax.plot(x_fit, y_fit, linestyle='--', linewidth=3, color=plt.gca().lines[-1].get_color(), label=parameter)

    ax.set_ylabel(r'$\left[ \eta \right]$', fontsize=20)
    ax.set_xlabel(r'$\Phi$', fontsize=20)
    ax.tick_params(which='both', labelsize=15, width=2, length=8, direction='in')
    ax.xaxis.set_major_locator(MaxNLocator(5))
    ax.xaxis.set_minor_locator(MaxNLocator(10))
    ax.yaxis.set_major_locator(MaxNLocator(5)) 
    ax.yaxis.set_minor_locator(MaxNLocator(10))
    #ax.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
    ax.set_title(plotTitle, fontsize=20)

    fig.tight_layout()
    plt.savefig('Pictures/DF_vs_IV_{}.png'.format(plotName), dpi=300, bbox_inches='tight')

if __name__ == "__main__":
    phi_range = [2.3994, 2.9993, 3.4492, 3.8991, 4.9488, 5.9986]
    Ca_range = [0.01, 0.02, 0.03, 0.06, 0.07, 0.08, 0.09, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2]

    '''
    # fix phi
    for phi in phi_range:
        ps = []
        for Ca in Ca_range:
            ps.append([phi, Ca])
        makeDFvsIVplot(ps, 1, '{} = {:.2}%'.format(r'$\phi$', phi), 'phi_{}'.format(phi))

    
    # fix Ca
    for Ca in Ca_range:
        ps = []
        for phi in phi_range:
            ps.append([phi, Ca])
        makeDFvsIVplot(ps, 1, 'Ca = {}'.format(Ca), 'Ca_{}'.format(Ca))
    
    '''
    ps = []
    for phi in [4.9488, 5.9986]:
        for Ca in [0.03, 0.06]:
            ps.append([phi, Ca])
    makeDFvsIVplot(ps, 1, '', 'manuscript')
    

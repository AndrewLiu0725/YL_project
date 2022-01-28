# ===============================================================================
# Copyright 2021 An-Jun Liu
# Last Modified Date: 01/23/2022
# ===============================================================================
import matplotlib.pyplot as plt
import numpy as np 
from RBC_Utilities import calcDoubletFraction, getIntrinsicViscosity, getSuspensionParameterSets

def makeDFvsIVplot(ps, plotType):
    
    plt.rcParams['font.family'] = 'DeJavu Serif'
    plt.rcParams['font.serif'] = ['Times New Roman']
    
    fig, ax = plt.subplots(figsize = (8, 6))
    [_, parameter_set] = getSuspensionParameterSets()
    st = 5

    for [phi, Ca] in ps:
        if plotType == 0:
            df_total, iv_total = [], []
        else:
            data = {}           

        for ensemble_id in parameter_set[phi][Ca]:
            try:
                [df, et, _] = calcDoubletFraction(phi, Ca, 1, 1, 0, ensemble_id, 1)
                iv = getIntrinsicViscosity(phi, Ca, None, ensemble_id, 1)

                if plotType == 0:
                    ax.scatter(df[0][st:et], iv[st:et], color='b', facecolor='none')
                    df_total += list(df[0][st:et])
                    iv_total += list(iv[st:et])
                else:
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
            ax.plot(df_total, iv_total, marker='s', markerfacecolor='None', linestyle='None', label='{}={:.2}%\nCa={}'.format(r'$\phi$', phi, Ca))

        model = np.polyfit(df_total, iv_total, 1)
        predict = np.poly1d(model)
        x_fit = np.linspace(min(df_total), max(df_total), 100)
        y_fit = predict(x_fit)
        ax.plot(x_fit, y_fit, linestyle='--', linewidth=3, color=plt.gca().lines[-1].get_color(), label='{}={:.2}%\nCa={}'.format(r'$\phi$', phi, Ca))

    ax.set_ylabel(r'$\left[ \eta \right]$', fontsize = 15)
    ax.set_xlabel(r'$\Phi$', fontsize = 15)
    ax.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
    #ax.set_title("phi = {}, Ca = {}".format(phi, Ca))
    #textstr = 'y={:.3}x+{:.3}'.format(model[0], model[1])
    #ax.text(0.05, 0.95, textstr, fontsize=15, transform=ax.transAxes)
    #plt.savefig('Pictures/DF_vs_IV_phi_{}_Ca_{}.png'.format(phi, Ca), dpi=200)
    fig.tight_layout()
    plt.savefig('Pictures/DF_vs_IV.png', dpi=300)

ps = []
for phi in [4.9488, 5.9986]:
    for Ca in [0.03, 0.06]:
        ps.append([phi, Ca])
makeDFvsIVplot(ps, 1)

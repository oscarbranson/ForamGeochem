import numpy as np
import matplotlib.pyplot as plt

from scipy import stats
from scipy.optimize import curve_fit
import uncertainties as un
from uncertainties.unumpy import nominal_values as nom
from uncertainties.unumpy import std_devs as err

from cbsyst import Csys

from .helpers import isolate_constant_conditions
from .plot import spreadm, angle_of_line

# Code for making supplementary Figures
#######################################

def Figure_SX(raw, dat, mdict, ldict):
    fig, ax = plt.subplots(1, 1)

    vlim = dat.loc[:, ('csys_mid', 'DIC')].min(), dat.loc[:, ('csys_mid', 'DIC')].max()

    for who in mdict.keys():
        ind = raw.Measured.who == who
        if who not in ['This Study', 'Haynes (subm)', 'Russell', 'Spero']:
            s = 25
            z = -1
            c = 'C3'
        else:
            s = 35
            z = 1
            c = 'C0'
        
        ma = ax.scatter(raw.loc[ind, ('csys_mid', 'CO3')], 
                        raw.loc[ind, ('csys_mid', 'pHtot')], 
                        marker=mdict[who], label=ldict[who],
                        color=c, alpha=0.75,
    #                     color=dat.loc[ind, ('csys_mid', 'DIC')],
                        vmin=vlim[0], vmax=vlim[1], cmap=plt.cm.Blues, 
                        edgecolor='k', lw=0.5, s=s, zorder=z)

    ax.set_xlabel('$[CO_3^{2-}]\ (\mu mol\ kg^{-1})$')
    ax.set_ylabel('$pH_{Total}$')
    # fig.colorbar(ma, label='[DIC]')

    ax.legend(loc='upper left', fontsize=8)

    ax.set_ylim(ax.get_ylim())
    ax.set_xlim(0, 550)

    fig.tight_layout()

    for DIC in [1000, 2000, 4000, 8000]:
        cs = Csys(np.linspace(*ax.get_ylim(), 50), DIC)
        line = ax.plot(cs.CO3, cs.pHtot, ls='dashed', color=(0,0,0,0.4), zorder=-1, lw=1)
        
        ind = (line[0].get_xdata() < ax.get_xlim()[1]) & (line[0].get_ydata() < ax.get_ylim()[1])
        
        x = line[0].get_xdata()[ind][-2:]
        y = line[0].get_ydata()[ind][-2:]
        angle = angle_of_line(x, y, ax)
        
        ax.text(x[0], y[0] + 0.01, 'DIC: {:.0f}'.format(DIC), rotation=angle, ha='right', alpha=0.4)
        
    return fig, ax


# Code for making figures in the Python notebook supplement.
############################################################

def parameter_space(dat):
    fig, ax = plt.subplots(1,1)

    x_thresh = [19, 21, 23, 25.5, 27]
    y_thresh = [1.8, 4, 8, 12]

    xt_last = 0
    yt_last = 0

    vmin = dat.loc[:, ('csys_mid', 'DIC')].min()
    vmax = dat.loc[:, ('csys_mid', 'DIC')].max()

    for xt in x_thresh:
        for yt in y_thresh:
            ind = ((dat.loc[:, ('Measured', 'Temp')] >= xt_last) & (dat.loc[:, ('Measured', 'Temp')] <= xt) &
                (dat.loc[:, ('Measured', 'Mg/Casw')] >= yt_last) & (dat.loc[:, ('Measured', 'Mg/Casw')] <= yt))
            if sum(ind) > 0:
                x, y, xm, ym = spreadm(dat.loc[ind, ('Measured', 'Temp')].astype(float).values,
                                       dat.loc[ind, ('Measured', 'Mg/Casw')].astype(float).values,
                                       x_tol=0.15, y_tol=0.17, offset_mult=0.2)

                if xm is not None:
                    for xi, yi in zip(x,y):
                        c = (.6, .6, .6)
                        ax.plot([xm, xi], [ym, yi], lw=0.5, color=c, zorder=-1)
                        ax.scatter(xm, ym, color='w', edgecolor=c, lw=0.5, s=10)

                ma = ax.scatter(x, y, c=sorted(dat.loc[ind, ('csys_mid', 'DIC')]), vmin=vmin, vmax=vmax, cmap=plt.cm.Blues, s=20, lw=0.5, edgecolor='k')       
            
            yt_last = yt
        xt_last = xt

    ax.set_xlabel('Temperature ($^{\circ}C$)')
    ax.set_ylabel('Mg/Ca$_{SW}$ (mol/mol)')
    
    fig.colorbar(ma, label='[DIC] ($\mu M$)')

    return fig, ax

def add_vs_mult(dat, crossplots):
    fig, axs = plt.subplots(2, 2)

    sub = isolate_constant_conditions(dat, Temp=22, tolerance=0.08)
    vmin, vmax = sub.loc[:, ('Measured', 'D_Mg')].min() , sub.loc[:, ('Measured', 'D_Mg')].max() * 1.2
    cmap = plt.cm.Blues_r

    ax = axs[0,0]

    ax.set_title('$D_{Mg} = C_1\ DIC + C_2\ [Ca] + D$', loc='left', fontsize=8)
    cm = ax.scatter(sub.loc[:, ('Measured', '[Ca]sw')], sub.loc[:, ('csys_mid', 'DIC')], c=sub.loc[:, ('Measured', 'D_Mg')] ,
                    edgecolors='k', lw=0.5, vmin=vmin, vmax=vmax, cmap=cmap)

    ax.set_xlabel('[Ca]sw')
    ax.set_ylabel('DIC')

    # fit an additive relationship
    def ca_dic_fn(x, C1, C2, D):
        Ca, DIC = x
        return  DIC * C1 + C2 * Ca + D

    p, cov = curve_fit(ca_dic_fn, (sub.loc[:, ('Measured', '[Ca]sw')], sub.loc[:, ('csys_mid', 'DIC')]), sub.loc[:, ('Measured', 'D_Mg')] )

    x = np.linspace(*ax.get_xlim(), 100)
    y = np.linspace(*ax.get_ylim(), 100)
    X, Y = np.meshgrid(x, y)

    ax.pcolormesh(X, Y, ca_dic_fn((X, Y), *p), zorder=-1, vmin=vmin, vmax=vmax, cmap=cmap)

    rax = axs[1, 0]
    rax.set_ylabel('N')
    resid = sub.loc[:, ('Measured', 'D_Mg')] - ca_dic_fn((sub.loc[:, ('Measured', '[Ca]sw')], sub.loc[:, ('csys_mid', 'DIC')]), *p)
    rax.hist(resid, bins=20)

    ax = axs[0,1]
    ax.set_title('$D_{Mg} = C_1\ DIC\ [Ca] + D$', loc='left', fontsize=8)
    cm = ax.scatter(sub.loc[:, ('Measured', '[Ca]sw')], sub.loc[:, ('csys_mid', 'DIC')], c=sub.loc[:, ('Measured', 'D_Mg')] ,
                    edgecolors='k', lw=0.5, vmin=vmin, vmax=vmax, cmap=cmap)
    ax.set_xlabel('[Ca]sw')

    # fit the relationship
    def ca_dic_fn(x, C1, D):
        Ca, DIC = x
        return  DIC * C1 * Ca + D

    p, cov = curve_fit(ca_dic_fn, (sub.loc[:, ('Measured', '[Ca]sw')], sub.loc[:, ('csys_mid', 'DIC')]), sub.loc[:, ('Measured', 'D_Mg')] )

    x = np.linspace(*ax.get_xlim(), 100)
    y = np.linspace(*ax.get_ylim(), 100)
    X, Y = np.meshgrid(x, y)

    ax.pcolormesh(X, Y, ca_dic_fn((X, Y), *p), zorder=-1, vmin=vmin, vmax=vmax, cmap=cmap)

    rax = axs[1, 1]

    resid = sub.loc[:, ('Measured', 'D_Mg')] - ca_dic_fn((sub.loc[:, ('Measured', '[Ca]sw')], sub.loc[:, ('csys_mid', 'DIC')]), *p)
    rax.hist(resid, bins=20)

    fig.tight_layout(rect=[0, 0, .85, 1])

    pos = ax.get_position()
    cax = fig.add_axes([.855, pos.y0, .02, pos.height])

    fig.colorbar(cm, cax=cax, label='$D_{Mg}$')


    for rax in axs[1, :]:
        rax.set_ylim(0, 11)
        rax.set_xlabel('Residual')
    
    crossplots['Ca_DIC'] = (X, Y, ca_dic_fn((X, Y), *p), sub.loc[:, ('Measured', '[Ca]sw')], sub.loc[:, ('csys_mid', 'DIC')], sub.loc[:, ('Measured', 'D_Mg')])

    return fig, axs, crossplots

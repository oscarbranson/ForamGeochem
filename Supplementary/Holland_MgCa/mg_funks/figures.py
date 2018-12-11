import numpy as np
from pandas import IndexSlice as idx
import matplotlib as mpl
import matplotlib.pyplot as plt

from matplotlib import ticker

from scipy import stats
from scipy.optimize import curve_fit
import uncertainties as un
from uncertainties.unumpy import nominal_values as nom
from uncertainties.unumpy import std_devs as err

from cbsyst import Csys

from .helpers import isolate_constant_conditions
from .plot import emphasise_subset, subset_trendline, contour_labels
from .model import T_fn

# Code for making figures in paper
##################################

def fig1(dat, rus, mdict, ldict):
    fig, ax = plt.subplots(1, 1)

    # plot russell
    russub = isolate_constant_conditions(rus, Mg=50, Ca=10.2, DIC=2000, pH=8.05)
    # ax.scatter(russub.loc[:, ('Measured', 'Temp')], russub.loc[:, ('Measured', 'Mg/Caf')], 
    #         marker='D', edgecolor='orange', facecolor='w', zorder=6, label='Russell et al (2004)')

    for who in dat.Measured.who.unique():
        whosub = dat.loc[dat.loc[:, ('Measured', 'who')] == who]
        subamb = isolate_constant_conditions(whosub, Mg=50, Ca=10.2, MgCa=5.1, DIC=2000, pH=8.05)
        subexp = whosub.drop(subamb.index)
        
        cation_exp = isolate_constant_conditions(subexp, DIC=2000, pH=8.05)
        anion_exp = isolate_constant_conditions(subexp, Mg=50, Ca=10.2, MgCa=5.1)
        both_var = subexp.drop(np.concatenate([cation_exp.index, anion_exp.index]))
    
        ax.errorbar(subamb.loc[:, ('Measured', 'Temp')],
                    subamb.loc[:, ('Measured', 'Mg/Caf')],
                    yerr=subamb.loc[:, ('Measured', 'Mg/Caf 2se')],
                    lw=0, elinewidth=1, marker=mdict[who], color='k', label=ldict[who], markersize=5)
        
        ax.errorbar(cation_exp.loc[:, ('Measured', 'Temp')],
                    cation_exp.loc[:, ('Measured', 'Mg/Caf')],
                    yerr=cation_exp.loc[:, ('Measured', 'Mg/Caf 2se')],
                    lw=0, elinewidth=1, marker=mdict[who], color='r', alpha=.6, markersize=5, label="_", zorder=2)
        
        ax.errorbar(anion_exp.loc[:, ('Measured', 'Temp')],
                    anion_exp.loc[:, ('Measured', 'Mg/Caf')],
                    yerr=anion_exp.loc[:, ('Measured', 'Mg/Caf 2se')],
                    lw=0, elinewidth=1, marker=mdict[who], color='b', alpha=.6, markersize=5, label="_", zorder=2)
        
        ax.errorbar(both_var.loc[:, ('Measured', 'Temp')],
                    both_var.loc[:, ('Measured', 'Mg/Caf')],
                    yerr=both_var.loc[:, ('Measured', 'Mg/Caf 2se')],
                    lw=0, elinewidth=1, marker=mdict[who], color='purple', alpha=.6, markersize=5, label="_", zorder=2)

    ax.scatter([], [], color='k', alpha = 0.3, marker='', label='Variable T$^{\circ}C$ only')
    ax.scatter([], [], color='r', alpha=0.3, marker='', label='Variable [Mg] and [Ca]')
    ax.scatter([], [], color='b', alpha=0.3, marker='', label='Variable pH and [DIC]')
    ax.scatter([], [], color=(.6, 0, .6), alpha=0.3, marker='', label='Multiple Variables')

    handles, labels = ax.get_legend_handles_labels()

    leg = ax.legend(handles[-3:] + handles[:-3], labels[-3:] + labels[:-3], loc='upper left', fontsize=8)

    for text in leg.get_texts():
        if '[Mg]' in text.get_text():
            text.set_color((1, 0, 0, 0.6))
        if '[DIC]' in text.get_text():
            text.set_color((0, 0, 1, 0.6))
        if 'Multiple' in text.get_text():
            text.set_color((.6, 0, .6, 0.6))

    def expfn(x, A, B):
        return B * np.exp(x * A)

    def expfn_err(x, A, B):
        return B * un.unumpy.exp(x * A)

    # whosub = dat.loc[dat.loc[:, ('Measured', 'who')] == 'This Study']
    ambient = isolate_constant_conditions(dat, Mg=50, Ca=10.2, MgCa=5.1, DIC=2000, pH=8.05, pH_tolerance=0.1)

    # draw ambient line
    p, cov = curve_fit(expfn, ambient.loc[:, ('Measured', 'Temp')], ambient.loc[:, ('Measured', 'Mg/Caf')])
    p_err = un.correlated_values(p, cov)

    # plot calibration and error envelope
    ax.set_xlim(ax.get_xlim())
    tn = np.linspace(*ax.get_xlim(), 100)
    pred = expfn_err(tn, *p_err)

    CI95_m = stats.t.interval(0.95, df=ambient.shape[0] - 1)[1]

    ax.plot(tn, nom(pred), zorder=-1, ls='dashed', lw=1, color='k')
    ax.fill_between(tn, 
                    nom(pred) - CI95_m * err(pred),
                    nom(pred) + CI95_m * err(pred),
                    zorder=-2, alpha=0.2, color='k', lw=0) 

    # Russell Line (ambient only)
    rp, rcov = curve_fit(expfn, russub.Measured.Temp, russub.loc[:, ('Measured', 'Mg/Caf')], 
                        sigma=russub.loc[:, ('Measured', 'Mg/Caf 2se')])
    rp_err = un.correlated_values(rp, rcov)

    rCI95_m = stats.t.interval(0.95, df=russub.shape[0] - 1)[1]

    rust = np.linspace(*ax.get_xlim())
    rusmgca = expfn_err(rust, *rp_err)
    ax.plot(rust, nom(rusmgca), ls='dashed', zorder=5, color='orange', lw=1)
    ax.fill_between(rust, 
                    nom(rusmgca) - rCI95_m * err(rusmgca),
                    nom(rusmgca) + rCI95_m * err(rusmgca), 
                    zorder=-2, alpha=0.3, color='orange', lw=0)


    ax.set_xlabel('Temperature ($^{\circ}C$)')
    ax.set_ylabel('$Mg/Ca_{\it{O. universa}}\ (mmol\ mol^{-1}$)')

    return fig, ax

def fig2(dat, mdict, ldict):

    labels = 'ABCDEFGHIJ'

    fig, axs = plt.subplots(2, 3, sharey=True)

    first = True

    # draw basic plots
    for k, m in mdict.items():
        ind = dat.Measured.who == k
        sub = dat.loc[ind, :]
        
        # Variable Temperature
        ax = axs[0,0]
        ax.scatter(sub.loc[:, ('Measured', 'Temp')], sub.loc[:, ('Measured', 'Mg/Caf')],
                s=sub.loc[:, ('Measured', 'numberforams')], marker=m, color=(0,0,0,0.1), lw=0)
        ax.set_xlim(ax.get_xlim())
        ax.set_xlabel('Temperature ($^{\circ}C$)')
        
        # Variable DIC
        ax = axs[0,1]
        ax.scatter(sub.loc[:, ('csys_mid', 'DIC')], sub.loc[:, ('Measured', 'Mg/Caf')],
                s=sub.loc[:, ('Measured', 'numberforams')], marker=m, color=(0,0,0,0.2), lw=0)
        ax.set_xlim(ax.get_xlim())
        ax.set_xlabel('DIC ($\mu mol\ kg^{-1}$)')

        # Variable pH
        ax = axs[0,2]
        ax.scatter(sub.loc[:, ('csys_mid', 'pHtot')], sub.loc[:, ('Measured', 'Mg/Caf')],
                s=sub.loc[:, ('Measured', 'numberforams')], marker=m, color=(0,0,0,0.2), lw=0)
        ax.set_xlim(ax.get_xlim())
        ax.set_xlabel('pH (Total)')

        # Mg/Ca variation 
        # Variable Mg
        ax = axs[1,0]
        ax.scatter(sub.loc[:, ('Measured', '[Mg]sw')], sub.loc[:, ('Measured', 'Mg/Caf')],
                s=sub.loc[:, ('Measured', 'numberforams')], marker=m, color=(0,0,0,0.2), lw=0)
        ax.set_xlim(ax.get_xlim())
        ax.set_xlabel('$[Mg]_{SW}\ (mmol\ kg^{-1})$')

        # Variable Ca
        ax = axs[1,1]
        ax.scatter(sub.loc[:, ('Measured', '[Ca]sw')], sub.loc[:, ('Measured', 'Mg/Caf')],
                s=sub.loc[:, ('Measured', 'numberforams')], marker=m, color=(0,0,0,0.2), lw=0)
        ax.set_xlim(ax.get_xlim())
        ax.set_xlabel('$[Ca]_{SW}\ (mmol\ kg^{-1})$')

        # Variable Mg/Ca
        ax = axs[1,2]
        ax.scatter(sub.loc[:, ('Measured', 'Mg/Casw')], sub.loc[:, ('Measured', 'Mg/Caf')],
                s=sub.loc[:, ('Measured', 'numberforams')], marker=m, color=(0,0,0,0.2), lw=0)
        ax.set_xlim(ax.get_xlim())
        ax.set_xlabel('$Mg/Ca_{SW}$')

    for ax in axs[:,0]:
        ax.set_ylabel('$Mg/Ca_{\it{O. universa}}$ (mmol $mol^{-1}$)')

    fig.tight_layout()

    first = True
    # draw subsets and trend lines
    for k, m in mdict.items():
        ind = dat.Measured.who == k
        sub = dat.loc[ind, :]

        # Variable Temperature   
        ax = axs[0,0]

        # ambient bros
        ssub = isolate_constant_conditions(sub, Ca=10, Mg=50, MgCa=5, DIC=2000, pH=8.1)
        emphasise_subset(ax, ssub, ('Measured', 'Temp'), color='w', m=m)
        if first:
            dsub = isolate_constant_conditions(dat, Ca=10, Mg=50, MgCa=5, DIC=2000, pH=8.1) 
            subset_trendline(ax, dsub, ('Measured', 'Temp'), color='k', line_order=1, 
                            label='Ambient', label_x=22.5, fontsize=6)
            
        # low Mg/Ca
        ssub = isolate_constant_conditions(sub, MgCa=1.3, DIC=2000, tolerance=0.2)
        emphasise_subset(ax, ssub, ('Measured', 'Temp'), color='C0', m=m)
        if first:
            dsub = isolate_constant_conditions(dat, MgCa=1.3, DIC=2000, tolerance=0.2) 
            subset_trendline(ax, dsub, ('Measured', 'Temp'), color='C0', line_order=1, 
                            label='$Mg/Ca_{SW} \\approx 1.3$', label_x=22.5, fontsize=6)

        # high Mg/Ca
        ssub = isolate_constant_conditions(sub, MgCa=10, DIC=2000)
        emphasise_subset(ax, ssub, ('Measured', 'Temp'), color='C3', m=m)
        if first:
            dsub = isolate_constant_conditions(dat, MgCa=10, DIC=2000) 
            subset_trendline(ax, dsub, ('Measured', 'Temp'), color='C3', line_order=1, 
                            label='$Mg/Ca_{SW} \\approx 10$', label_x=22.5, fontsize=6)

        # Variable DIC
        ax = axs[0,1]

        # ambient bros, low-pH
        ssub = isolate_constant_conditions(sub, Temp=22, Ca=10, Mg=50, MgCa=5)
        emphasise_subset(ax, ssub, ('csys_mid', 'DIC'), color=(.2,.2,.2), m='x', s=15)
        if first:
            dsub = isolate_constant_conditions(dat, Temp=22, Ca=10, Mg=50, MgCa=5) 
            subset_trendline(ax, dsub, ('csys_mid', 'DIC'), color='k', line_order=1, 
                            label='Ambient\n(low pH)', label_x=5000, fontsize=6, label_v_offset_frac=0, lsargs={'ls': 'dashed'})

        # ambient bros
        ssub = isolate_constant_conditions(sub, Temp=22, Ca=10, Mg=50, MgCa=5, pH=8.1)
        emphasise_subset(ax, ssub, ('csys_mid', 'DIC'), color='w', m=m)

        # low Mg/Ca
        ssub = isolate_constant_conditions(sub, MgCa=1.3, Temp=22, tolerance=0.2)
        emphasise_subset(ax, ssub, ('csys_mid', 'DIC'), color='C0', m=m)
        if first:
            dsub = isolate_constant_conditions(dat, MgCa=1.3, Temp=22, tolerance=0.2) 
            subset_trendline(ax, dsub, ('csys_mid', 'DIC'), color='C0', line_order=1, 
                            label=None, label_x=5000, fontsize=6)
            
        # high Mg/Ca
        ssub = isolate_constant_conditions(sub, MgCa=10, Temp=22)
        emphasise_subset(ax, ssub, ('csys_mid', 'DIC'), color='C3', m=m)
        if first:
            dsub = isolate_constant_conditions(dat, MgCa=10, Temp=22) 
            subset_trendline(ax, dsub, ('csys_mid', 'DIC'), color='C3', line_order=1, 
                            label=None, label_x=5000, fontsize=6)

        # Variable pH
        ax = axs[0,2]

        # ambient bros
        ssub = isolate_constant_conditions(sub, Temp=22, MgCa=5.1, DIC=2000)
        emphasise_subset(ax, ssub, ('csys_mid', 'pHtot'), color='w', m=m)
        if first:
            dsub = isolate_constant_conditions(dat, Temp=22, MgCa=5.1, DIC=2000) 
            subset_trendline(ax, dsub, ('csys_mid', 'pHtot'), color='k', line_order=1, 
                            label=None, label_x=8.2, fontsize=6)

        # low Mg/Ca
        ssub = isolate_constant_conditions(sub, MgCa=1.3, Temp=22, DIC=2000, tolerance=0.2)
        emphasise_subset(ax, ssub, ('csys_mid', 'pHtot'), color='C0', m=m)
        if first:
            dsub = isolate_constant_conditions(dat, MgCa=1.3, Temp=22, DIC=2000, tolerance=0.2) 
            subset_trendline(ax, dsub, ('csys_mid', 'pHtot'), color='C0', line_order=1, 
                            label=None, label_x=8.2, fontsize=6)

        # low Mg/Ca
        ssub = isolate_constant_conditions(sub, MgCa=10, Temp=22, DIC=2000, tolerance=0.1)
        emphasise_subset(ax, ssub, ('csys_mid', 'pHtot'), color='C3', m=m)
        if first:
            dsub = isolate_constant_conditions(dat, MgCa=10, Temp=22, DIC=2000, tolerance=0.1) 
            subset_trendline(ax, dsub, ('csys_mid', 'pHtot'), color='C3', line_order=1, 
                            label=None, label_x=8.2, fontsize=6)

        # Mg/Ca variation 
        # Variable Mg
        ax = axs[1,0]

        # ambient bros
        ssub = isolate_constant_conditions(sub, Temp=22, Ca=10, Mg=None, MgCa=None, DIC=2000, pH=8.1)
        emphasise_subset(ax, ssub, ('Measured', '[Mg]sw'), color='w', m=m)
        if first:
            dsub = isolate_constant_conditions(dat, Temp=22, Ca=10, Mg=None, MgCa=None, DIC=2000, pH=8.1) 
            subset_trendline(ax, dsub, ('Measured', '[Mg]sw'), color='k', line_order=1, 
                            label='$[Ca] \\approx 10$', label_x=65, fontsize=6)

        # Variable Ca
        ax = axs[1,1]

        # ambient bros
        ssub = isolate_constant_conditions(sub,  Temp=22, Ca=None, Mg=50, MgCa=None, DIC=None, pH=None)
        emphasise_subset(ax, ssub, ('Measured', '[Ca]sw'), color='w', m=m)
        if first:
            dsub = isolate_constant_conditions(dat,  Temp=22, Ca=None, Mg=50, MgCa=None, DIC=None, pH=None) 
            subset_trendline(ax, dsub, ('Measured', '[Ca]sw'), color='k', line_order=1, 
                            label='$[Mg] \\approx 50$', label_x=12, fontsize=6)
        # low [Mg]
        ssub = isolate_constant_conditions(sub, Temp=22, Ca=None, Mg=25, MgCa=None, DIC=None, pH=None)
        emphasise_subset(ax, ssub, ('Measured', '[Ca]sw'), color='C0', m=m)
        if first:
            dsub = isolate_constant_conditions(dat, Temp=22, Ca=None, Mg=25, MgCa=None, DIC=None, pH=None) 
            subset_trendline(ax, dsub, ('Measured', '[Ca]sw'), color='C0', line_order=1, 
                            label='$[Mg] \\approx 25$', label_x=21, fontsize=6, label_v_offset_frac=-0.05)

        # high [Mg]
        ssub = isolate_constant_conditions(sub, Temp=22, Ca=None, Mg=100, MgCa=None, DIC=None, pH=None)
        emphasise_subset(ax, ssub, ('Measured', '[Ca]sw'), color='C3', m=m)
        ax.text(12, 9.5, '$[Mg] \\approx 100$', color='C3', ha='left', fontsize=6, alpha=0.3)


        # Variable Mg/Ca
        ax = axs[1,2]
        # ambient bros
        ssub = isolate_constant_conditions(sub,  Temp=22, Ca=None, Mg=None, MgCa=None, DIC=2000, pH=8.1)
        emphasise_subset(ax, ssub, ('Measured', 'Mg/Casw'), 'w', m=m)
        if first:
            dsub = isolate_constant_conditions(dat,  Temp=22, Ca=None, Mg=None, MgCa=None, DIC=2000, pH=8.1) 
            subset_trendline(ax, dsub, ('Measured', 'Mg/Casw'), color='k', line_order=1, 
                            label=None, label_x=7.7, fontsize=6, label_v_offset_frac=0.05)

        first = False

    # size scale bar
    ax = axs[-1,-1]
    for s in [dat.loc[:, ('Measured', 'numberforams')].min(),
            dat.loc[:, ('Measured', 'numberforams')].median(),
            dat.loc[:, ('Measured', 'numberforams')].max()]:

        ax.scatter([], [], s=s, color='w', edgecolor='k', marker='o', lw=0.5, label='{:.0f}'.format(s))

    hans, labs = ax.get_legend_handles_labels()
    leg = ax.legend(hans[-3:], labs[-3:], loc='lower right', fontsize=8)

    for i, ax in enumerate(axs.flat):
        ax.text(.03, .97, labels[i], ha='left', va='top', transform=ax.transAxes, 
                fontsize=12, weight='bold', color=(.3,.3,.3), zorder=2)

    return fig, axs

def fig3(dat, mdict):
    labels = 'ABCDEFGHIJ'

    fig, axs = plt.subplots(2, 3, sharey=True)

    for k, m in mdict.items():
        ind = dat.Measured.who == k
        sub = dat.loc[ind, :]
        # Variable Temperature
        ax = axs[0,0]
        ax.scatter(sub.loc[:, ('Measured', 'Temp')], sub.loc[:, ('Measured', 'D_Mg')],
                s=sub.loc[:, ('Measured', 'numberforams')], marker=m, color=(0,0,0,0.1), lw=0)
        ax.set_xlim(ax.get_xlim())
        ax.set_xlabel('Temperature ($^{\circ}C$)')
        
        # Variable DIC
        ax = axs[0,1]
        ax.scatter(sub.loc[:, ('csys_mid', 'DIC')], sub.loc[:, ('Measured', 'D_Mg')],
                s=sub.loc[:, ('Measured', 'numberforams')], marker=m, color=(0,0,0,0.2), lw=0)
        ax.set_xlim(ax.get_xlim())
        ax.set_xlabel('DIC ($\mu mol\ kg^{-1}$)')

        # Variable pH
        ax = axs[0,2]
        ax.scatter(sub.loc[:, ('csys_mid', 'pHtot')], sub.loc[:, ('Measured', 'D_Mg')],
                s=sub.loc[:, ('Measured', 'numberforams')], marker=m, color=(0,0,0,0.2), lw=0)
        ax.set_xlim(ax.get_xlim())
        ax.set_xlabel('pH (Total)')

        # Variable Mg
        ax = axs[1,0]
        ax.scatter(sub.loc[:, ('Measured', '[Mg]sw')], sub.loc[:, ('Measured', 'D_Mg')],
                s=sub.loc[:, ('Measured', 'numberforams')], marker=m, color=(0,0,0,0.2), lw=0)
        ax.set_xlim(ax.get_xlim())
        ax.set_xlabel('$[Mg]_{SW}\ (mmol\ kg^{-1})$')

    for ax in axs[:,0]:
        ax.set_ylabel('$D_{Mg} \\times 1000$')

    fig.tight_layout()

    # Draw subsets and trendlines
    first = True
    for k, m in mdict.items():
        ind = dat.Measured.who == k
        sub = dat.loc[ind, :]

        # Variable Temperature
        ax = axs[0,0]
        # ambient
        ssub = isolate_constant_conditions(sub, Ca=10, DIC=2000, pH=8.1)
        emphasise_subset(ax, ssub, ('Measured', 'Temp'), color='w', yvar=('Measured', 'D_Mg'), m=m)
        if first:
            dsub = isolate_constant_conditions(dat, Ca=10, DIC=2000, pH=8.1)
            subset_trendline(ax, dsub, ('Measured', 'Temp'), 'k', 1, label='Ambient', label_x=22.5, 
                            yvar=('Measured', 'D_Mg'), fontsize=6)

        # high [Ca]
        ssub = isolate_constant_conditions(sub, Ca=20, DIC=2000, pH=None)
        emphasise_subset(ax, ssub, ('Measured', 'Temp'), color='C1', yvar=('Measured', 'D_Mg'), m=m)
        if first:
            dsub = isolate_constant_conditions(dat, Ca=20, DIC=2000, pH=None)
            subset_trendline(ax, dsub, ('Measured', 'Temp'), 'C1', 1, label='$[Ca]_{SW} \\approx 20$', label_x=22.5,
                            yvar=('Measured', 'D_Mg'), fontsize=6)

        # intermediate [Ca]
        ssub = isolate_constant_conditions(sub, Ca=15, DIC=2000, pH=None)
        emphasise_subset(ax, ssub, ('Measured', 'Temp'), color='C6', yvar=('Measured', 'D_Mg'), m=m)
        if first:
            dsub = isolate_constant_conditions(dat, Ca=15, DIC=2000, pH=None)
            subset_trendline(ax, dsub, ('Measured', 'Temp'), 'C6', 1, label='$[Ca]_{SW} \\approx 15$', label_x=22.5,
                            yvar=('Measured', 'D_Mg'), fontsize=6)        

        # Variable DIC
        ax = axs[0,1]
        # ambient bros, low-pH
        ssub = isolate_constant_conditions(sub, Temp=22, Ca=10, Mg=50, MgCa=5)
        emphasise_subset(ax, ssub, ('csys_mid', 'DIC'), color=(.2,.2,.2), yvar=('Measured', 'D_Mg'), m='x', s=15)
        if first:
            dsub = isolate_constant_conditions(dat, Temp=22, Ca=10, Mg=50, MgCa=5)
            subset_trendline(ax, dsub, ('csys_mid', 'DIC'), color='k', line_order=1, yvar=('Measured', 'D_Mg'), 
                            label='Ambient\n(low pH)', label_x=5000, fontsize=6, label_v_offset_frac=0, 
                                lsargs={'ls': 'dashed'})

        # ambient
        ssub = isolate_constant_conditions(sub, Ca=10, Temp=22, pH=8.1)
        emphasise_subset(ax, ssub, ('csys_mid', 'DIC'), color='w', yvar=('Measured', 'D_Mg'), m=m)

        # high [Ca]
        ssub = isolate_constant_conditions(sub, Ca=20, Temp=22, pH=None)
        emphasise_subset(ax, ssub, ('csys_mid', 'DIC'), color='C1', yvar=('Measured', 'D_Mg'), m=m)
        if first:
            dsub = isolate_constant_conditions(dat, Ca=20, Temp=22, pH=None)
            subset_trendline(ax, dsub, ('csys_mid', 'DIC'), 'C1', 1, label=None, label_x=22.5,
                            yvar=('Measured', 'D_Mg'), fontsize=6)

        # low [Ca]
        ssub = isolate_constant_conditions(sub, Ca=5, Temp=22, pH=None)
        emphasise_subset(ax, ssub, ('csys_mid', 'DIC'), color='C4', yvar=('Measured', 'D_Mg'), m=m)
        if first:
            dsub = isolate_constant_conditions(dat, Ca=5, Temp=22, pH=None)
            subset_trendline(ax, dsub, ('csys_mid', 'DIC'), 'C4', 1, label='$[Ca]_{SW} \\approx 5$', label_x=3000,
                            yvar=('Measured', 'D_Mg'), fontsize=6, label_v_offset_frac=-0.08)         

        # Variable pH
        ax = axs[0,2]

        # ambient bros
        ssub = isolate_constant_conditions(sub, Temp=22, Ca=10, DIC=2000)
        emphasise_subset(ax, ssub, ('csys_mid', 'pHtot'), color='w', yvar=('Measured', 'D_Mg'), m=m)
        if first:
            dsub = isolate_constant_conditions(dat, Temp=22, Ca=10, DIC=2000)
            subset_trendline(ax, dsub, ('csys_mid', 'pHtot'), 'k', 1,
                            yvar=('Measured', 'D_Mg'))

        # high [Ca]
        ssub = isolate_constant_conditions(sub, Ca=20, Temp=22, pH=None, DIC=2000)
        emphasise_subset(ax, ssub, ('csys_mid', 'pHtot'), color='C1', yvar=('Measured', 'D_Mg'), m=m)
        # if first:
        #     dsub = isolate_constant_conditions(dat, Ca=20, Temp=22, pH=None, DIC=2000)
        #     subset_trendline(ax, dsub, ('csys_mid', 'pHtot'), 'C1', 1, label=None, label_x=22.5,
        #                     yvar=('Measured', 'D_Mg'), fontsize=6)
        

        # low [Ca]
        ssub = isolate_constant_conditions(sub, Ca=5, Temp=22, pH=None, DIC=2000)
        emphasise_subset(ax, ssub, ('csys_mid', 'pHtot'), color='C4', yvar=('Measured', 'D_Mg'), m=m)
        if first:
            dsub = isolate_constant_conditions(dat, Ca=5, Temp=22, pH=None, DIC=2000)
            subset_trendline(ax, dsub, ('csys_mid', 'pHtot'), 'C4', 1, label=None, label_x=3000,
                            yvar=('Measured', 'D_Mg'), fontsize=6, label_v_offset_frac=-0.08) 

        # Variable Mg
        ax = axs[1,0]

        # ambient
        ssub = isolate_constant_conditions(sub, Temp=22, Ca=10, DIC=2000, pH=8.1)
        emphasise_subset(ax, ssub, ('Measured', '[Mg]sw'), color='w', yvar=('Measured', 'D_Mg'), m=m)
        if first:
            dsub = isolate_constant_conditions(dat, Temp=22, Ca=10, DIC=2000, pH=8.1)
            subset_trendline(ax, dsub, ('Measured', '[Mg]sw'), 'k', 1,
                             yvar=('Measured', 'D_Mg'), label='$[Ca] \\approx 10$', label_x=65, fontsize=6)

        # Variable Ca
        ax = axs[1,1]
        ax.scatter(sub.loc[:, ('Measured', '[Ca]sw')], sub.loc[:, ('Measured', 'D_Mg')],
                    s=sub.loc[:, ('Measured', 'numberforams')], marker=m, color=(0,0,0,0.2), lw=0)
        ax.set_xlim(ax.get_xlim())
        ax.set_xlabel('$[Ca]_{SW}\ (mmol\ kg^{-1})$')

        # ambient
        ssub = isolate_constant_conditions(sub, Temp=22, DIC=2000, pH=8.1)
        emphasise_subset(ax, ssub, ('Measured', '[Ca]sw'), color='w', yvar=('Measured', 'D_Mg'), m=m)
        if first:
            dsub = isolate_constant_conditions(dat, Temp=22, DIC=2000, pH=8.1)
            subset_trendline(ax, dsub, ('Measured', '[Ca]sw'), 'k', 1,
                            yvar=('Measured', 'D_Mg'))

        # Low [DIC]
        ssub = isolate_constant_conditions(sub, Temp=22, DIC=1000, pH=None)
        emphasise_subset(ax, ssub, ('Measured', '[Ca]sw'), color='C2', yvar=('Measured', 'D_Mg'), m=m)
        if first:
            dsub = isolate_constant_conditions(dat, Temp=22, DIC=1000, pH=None)
            subset_trendline(ax, dsub, ('Measured', '[Ca]sw'), 'C2', 1,
                            yvar=('Measured', 'D_Mg'), label='$DIC \\approx 1000$', label_x=21, fontsize=6,
                                label_v_offset_frac=-0.08)

        # High [DIC]
        ssub = isolate_constant_conditions(sub, Temp=22, DIC=4000, pH=None)
        emphasise_subset(ax, ssub, ('Measured', '[Ca]sw'), color='C0', yvar=('Measured', 'D_Mg'), m=m)
        if first:
            dsub = isolate_constant_conditions(dat, Temp=22, DIC=4000, pH=None)
            subset_trendline(ax, dsub, ('Measured', '[Ca]sw'), 'C0', 1,
                            yvar=('Measured', 'D_Mg'), label='$DIC \\approx 4000$', label_x=19, fontsize=6, 
                                label_v_offset_frac=0.08)

        # High T
        ssub = isolate_constant_conditions(sub, Temp=26, DIC=2000, pH=None)
        emphasise_subset(ax, ssub, ('Measured', '[Ca]sw'), color='C5', yvar=('Measured', 'D_Mg'), m=m)
        if first:
            dsub = isolate_constant_conditions(dat, Temp=26, DIC=2000, pH=None)
            subset_trendline(ax, dsub, ('Measured', '[Ca]sw'), 'C5', 1,
                            yvar=('Measured', 'D_Mg'), label='$Temp \\approx 26$', label_x=13,
                                fontsize=6)
        first = False

    axs[1, 2].set_visible(False)

    for i, ax in enumerate(axs.flat):
        ax.text(.03, .97, labels[i], ha='left', va='top', transform=ax.transAxes, 
                fontsize=12, weight='bold', color=(.3,.3,.3), zorder=2)

    return fig, axs

def fig4(dat, params, pred, pred_T, ifa_rsd, mdict, cis_upper, cis_lower, med_ifa):
    labels = 'ABCDEFGHIJ'

    fig, (mgax, tax) = plt.subplots(1, 2, figsize=[7.5, 3])

    opts = {'cmap': plt.cm.Blues,
            'lw': 0.5,
            'edgecolor': 'k'}

    eopts = {'lw': 0,
            'elinewidth': 1,
            'color': (0,0,0,0.3),
            'zorder': -1}

    for s in ['This Study', 'Spero', 'Russell']:
        ind = dat.Measured.who == s
        mgax.scatter(dat.loc[ind, ('Measured', 'Mg/Caf')], nom(pred[ind]), marker=mdict[s],
                    s=dat.loc[ind, ('Measured', 'numberforams')], c=dat.loc[ind, ('csys_mid', 'pHtot')], **opts)
        
        mgax.errorbar(dat.loc[ind, ('Measured', 'Mg/Caf')], nom(pred[ind]),
                    xerr=dat.loc[ind, ('Uncertainties', 'CI95')], yerr=1.96 * err(pred[ind]),
                    **eopts)
        
        cb = tax.scatter(dat.loc[ind, ('Measured', 'Temp')], nom(pred_T[ind]), marker=mdict[s],
                        s=dat.loc[ind, ('Measured', 'numberforams')], c=dat.loc[ind, ('csys_mid', 'pHtot')], **opts)
        tax.errorbar(dat.loc[ind, ('Measured', 'Temp')], nom(pred_T[ind]),
                    yerr=1.96 * err(pred_T[ind]),
                    **eopts)

    lim = np.array(mgax.get_xlim())
    mgax.set_xlim(lim)
    mgax.set_ylim(lim)
    mgax.plot(lim, lim, zorder=-2, ls='dashed', c=(0,0,0,0.5))

    for n in [4, 10, 30]:
        i = stats.t.interval(0.95, df=n-1)[-1]
        if n == 4:
            lab = '95% CI'
        else:
            lab = '_'
        mgax.fill_between(lim, 
                        lim + i * (lim * ifa_rsd) / n**0.5, 
                        lim - i * (lim * ifa_rsd) / n**0.5, 
                        alpha=0.1, color='k', lw=0, zorder=-1, label=lab)

    mgax.text(16.2, 9.7, 'N=4',  ha='right', rotation=28, fontsize=7, alpha=0.5)
    mgax.text(16.2, 13.2, 'N=10', ha='right', rotation=36, fontsize=7, alpha=0.5)
    mgax.text(16.2, 14.6, 'N=30', ha='right', rotation=38, fontsize=7, alpha=0.5)

    mgax.legend(loc=(0.02, 0.75), framealpha=0)


    lim = np.array((13, 30))
    tax.set_xlim(lim)
    tax.set_ylim(lim)
    tax.plot(lim, lim, zorder=-2, ls='dashed', c=(0,0,0,0.5))

    for upper_ifa, lower_ifa in zip(cis_upper, cis_lower):
        tax.fill_between(med_ifa, upper_ifa, lower_ifa, alpha=0.1, color='k', lw=0, zorder=-3, label='_')
        
    tax.text(30.0, 24, 'N=4',  ha='right', rotation=40, fontsize=7, alpha=0.5)
    tax.text(30.0, 27.5, 'N=10', ha='right', rotation=40, fontsize=7, alpha=0.5)
    tax.text(30.0, 28.5, 'N=30', ha='right', rotation=40, fontsize=7, alpha=0.5)    

    mgax.set_xlabel('Measured $Mg/Ca_{\it{O. universa}}$ (mmol $mol^{-1}$)')
    mgax.set_ylabel('Predicted $Mg/Ca_{\it{O. universa}}$ (mmol $mol^{-1}$)')

    tax.set_xlabel('Measured Temperature ($^{\circ}C$)')
    tax.set_ylabel('Predicted Tempterature ($^{\circ}C$)')

    label1 = '$Mg/Ca_{\it{O. universa}} = Mg/Ca_{SW}^A B\ e^{(C_1\ [Ca]_{SW} + C_2\ [DIC]_{SW} + D) T}$'

    label2 = ("A: {:.3f}\n".format(params[0]) +
            "B: {:.3f}\n".format(params[1]) +
            "C1: {:.3f}\n".format(params[2]) +
            "C2: {:.3f}\n".format(params[3]) +
            "D: {:.3f}\n".format(params[4]))
    fig.tight_layout(rect=[0,0,.9,1])

    mgax.text(.98, .02, label1, ha='right', va='bottom', transform=mgax.transAxes, fontsize=6)
    mgax.text(.98, .06, label2, ha='right', va='bottom', transform=mgax.transAxes, fontsize=6, color=(.3,.3,.3))

    # size scale bar
    for s in [4, 10, 30]:
        tax.scatter([], [], s=s, c='w', edgecolor='k', marker='o', lw=0.5, label='{:.0f}'.format(s))
    hans, labs = tax.get_legend_handles_labels()
    leg = tax.legend(hans[-3:], labs[-3:], loc='lower right', fontsize=8, framealpha=0.5)

    for i, ax in enumerate([mgax, tax]):
        ax.text(.03, .97, labels[i], ha='left', va='top', transform=ax.transAxes, 
                fontsize=14, weight='bold', color=(.3,.3,.3), zorder=2)

    cbax = fig.add_axes([tax.get_position().x1 + 0.025, tax.get_position().y0, .02, tax.get_position().height])
    fig.colorbar(cb, cax=cbax, label='$pH_{Total}$')

    # ticks
    mgax.xaxis.set_major_locator(ticker.MultipleLocator(5))
    mgax.yaxis.set_major_locator(ticker.MultipleLocator(5))
    mgax.xaxis.set_minor_locator(ticker.MultipleLocator(1))
    mgax.yaxis.set_minor_locator(ticker.MultipleLocator(1))
    
    tax.xaxis.set_major_locator(ticker.MultipleLocator(5))
    tax.yaxis.set_major_locator(ticker.MultipleLocator(5))
    tax.xaxis.set_minor_locator(ticker.MultipleLocator(1))
    tax.yaxis.set_minor_locator(ticker.MultipleLocator(1))

    return fig, ax

def fig5(rd, gray, rmdict, rrdict, rpred, rpar, tstyle):
    fig, axs = plt.subplots(2, 2, figsize=[6, 6], sharey=True, sharex=True)

    axs = axs.flatten()

    lim = (2, 9)

    labels = 'ABCDE'

    # draw DIC fit
    ax = axs[0]
    t = 'white'
    pt = 'white_DIC'
    d = rd.loc[rd.Type == t]

    for w in d.who.unique():
        cb = ax.scatter(d.loc[d.who == w, 'Mg/Caf'], nom(rpred[pt][d.who == w]), 
                        marker=rmdict[w], label=rrdict[w], s=15, edgecolor='k', lw=0.5, cmap=plt.cm.Blues,
                        c=d.loc[d.who == w, 'CO3'], vmin=75, vmax=530)

        ax.errorbar(d.loc[d.who == w, 'Mg/Caf'], nom(rpred[pt][d.who == w]), 
                    yerr=err(rpred[pt][d.who == w]),
                    xerr=d.loc[d.who == w, 'Mg/Caf 2se'], lw=0, elinewidth=1, color=(0,0,0,0.3), zorder=-1)


    for i, t in zip([1, 3], ['white', 'pink']):
        ax = axs[i]
        d = rd.loc[rd.Type == t]
        
        for w in d.who.unique():
            ax.scatter(d.loc[d.who == w, 'Mg/Caf'], nom(rpred[t][d.who == w]), 
                    marker=rmdict[w], label=rrdict[w], s=15, **tstyle[t])
            
            ax.errorbar(d.loc[d.who == w, 'Mg/Caf'], nom(rpred[t][d.who == w]), 
                        yerr=err(rpred[t][d.who == w]),
                        xerr=d.loc[d.who == w, 'Mg/Caf 2se'], lw=0, elinewidth=1, color=(0,0,0,0.3), zorder=-1)
        
        ax.legend(fontsize=7, loc='lower right', framealpha=0.2, borderpad=0.3)
        
        upe = rpar[t]
        upeph = rpar[t + '_pH']
        if t == 'white':
            parlab = ("A: {:.3f}\n".format(upe[0]) +
                    "B: {:.3f}\n".format(upe[1]) +
                    "$C_1$" + ": {:.1f}\n".format(upe[2]) +
                    "$C_2$" + ": {:.3f}\n".format(upe[3]) + 
                    "D: {:.3f}\n".format(upe[4]))
        else:
            parlab = ("P: {:.3f}\n".format(upe[0]) +
            "$C_1$" + ": {:.1f}\n".format(upe[1]) +
            "$C_2$" + ": {:.3f}\n".format(upe[2]) + 
            "D: {:.3f}\n".format(upe[3]))

        y = 0.76
        # if i == 3:
        #     y = 0.76
        # else:
        #     y = 0.84
        ax.text(.03, y, parlab, va='top', ha='left', fontsize=6, color=(.3, .3, .3), transform=ax.transAxes)
        
    # Gray data
    ax = axs[2]
    ax.scatter(gray.loc[:, 'G. ruber w Mg/Ca [mmol/mol]'], nom(rpred['gray']), 
            marker='.', label='Gray et al. (2018)', s=15, **tstyle['white'], alpha=0.6)
    ax.legend(fontsize=7, loc='lower right', framealpha=0.2, borderpad=0.3)

    upe = rpar['gray']
    upeph = rpar['gray_pH']
    parlab = ("A: {:.3f}\n".format(upe[0]) +
            "B: {:.3f}\n".format(upe[1]) +
            "$C_1$" + ": {:.1f}\n".format(upe[2]) +
            "$C_2$" + ": {:.3f}\n".format(upe[3]) + 
    #           "$C_2 (H)$" + ": {:.3f}$\\times 10^6$\n".format(upeph[3] * 1e-6) + 
            "D: {:.3f}\n".format(upe[4]))
    ax.text(.05, .84, parlab, va='top', ha='left', fontsize=6, color=(.3, .3, .3), transform=ax.transAxes)

    desc = ['culture - white', 'culture - white', 'sed. trap - white', 'culture - pink']
    # axis labels
    for i, ax in enumerate(axs):
        ax.set_xlim(lim)
        ax.set_ylim(lim[0], lim[-1] * 1.15)
        ax.plot(lim, lim, zorder=-1, ls='dashed', color=(0,0,0,0.4))
        
        if ax.is_last_row():
            ax.set_xlabel('Measured $Mg/Ca_{\it{G. ruber}}$ (mmol $mol^{-1}$)')
        if ax.is_first_col():
            ax.set_ylabel('Predicted $Mg/Ca_{\it{G. ruber}}$ (mmol $mol^{-1}$)')
        ax.text(.05, .95, labels[i], ha='left', va='top', transform=ax.transAxes, weight='bold',
                fontsize=14, color=(.3, .3, .3), zorder=2)
        ax.text(.125, .935, desc[i], ha='left', va='top', transform=ax.transAxes, weight='bold', style='italic',
                fontsize=8, color=(.6, .6, .6), zorder=2)

    ax = axs[1]
    eqnlab = '$Mg/Ca_{foram} = Mg/Ca_{SW}^A\ B\ e^{(C_1\ [Ca]_{SW}\ + C_2\ [CO_3]_{SW} + D) T}$\n'
    ax.text(.03, .86, eqnlab, va='top', ha='left', fontsize=7, color=(.15, .15, .15), transform=ax.transAxes)
    
    ax = axs[3]
    eqnlab = '$Mg/Ca_{foram} = P\ e^{(C_1\ [Ca]_{SW}\ + C_2\ [CO_3]_{SW} + D) T}$\n'
    ax.text(.03, .86, eqnlab, va='top', ha='left', fontsize=7, color=(.15, .15, .15), transform=ax.transAxes)

    fig.tight_layout(h_pad=0.01)

    for i, p in zip([1,3], ['white', 'pink']):
        upe = rpar[p]
        
    # CO3 colour scale
    pos = axs[0].get_position()
    cbax = fig.add_axes([pos.x1 - pos.width * 0.45, pos.y0 + pos.height * 0.15,
                        pos.width * 0.4, pos.height * 0.05])
    fig.colorbar(cb, cax=cbax, orientation='horizontal', ticks=[100, 300, 500], label='$[CO_3^{2-}]$')
    cbax.xaxis.set_label_position('top')

    # ticks
    for ax in axs:
        ax.xaxis.set_major_locator(ticker.MultipleLocator(2))
        ax.yaxis.set_major_locator(ticker.MultipleLocator(2))
        ax.xaxis.set_minor_locator(ticker.MultipleLocator(.5))
        ax.yaxis.set_minor_locator(ticker.MultipleLocator(.5))
    

    return fig, axs

def fig5_pH(rd, gray, rmdict, rrdict, rpred, rpar, tstyle):
    fig, axs = plt.subplots(2, 2, figsize=[6, 6], sharey=True, sharex=True)

    axs = axs.flatten()

    lim = (2, 9)
    labels = 'ABCDE'

    # draw DIC fit
    ax = axs[0]
    t = 'white'
    pt = 'white_DIC'
    d = rd.loc[rd.Type == t]

    for w in d.who.unique():
        cb = ax.scatter(d.loc[d.who == w, 'Mg/Caf'], nom(rpred[pt][d.who == w]), 
                        marker=rmdict[w], label=rrdict[w], s=15, edgecolor='k', lw=0.5, cmap=plt.cm.Blues,
                        c=d.loc[d.who == w, 'pHtot'], vmin=7.5, vmax=8.5)

        ax.errorbar(d.loc[d.who == w, 'Mg/Caf'], nom(rpred[pt][d.who == w]), 
                    yerr=err(rpred[pt][d.who == w]),
                    xerr=d.loc[d.who == w, 'Mg/Caf 2se'], lw=0, elinewidth=1, color=(0,0,0,0.3), zorder=-1)


    for i, t in zip([1, 3], ['white', 'pink']):
        ax = axs[i]
        d = rd.loc[rd.Type == t]
        
        for w in d.who.unique():
            ax.scatter(d.loc[d.who == w, 'Mg/Caf'], nom(rpred[t + '_pH'][d.who == w]), 
                    marker=rmdict[w], label=rrdict[w], s=15, **tstyle[t])
            
            ax.errorbar(d.loc[d.who == w, 'Mg/Caf'], nom(rpred[t + '_pH'][d.who == w]), 
                        yerr=err(rpred[t + '_pH'][d.who == w]),
                        xerr=d.loc[d.who == w, 'Mg/Caf 2se'], lw=0, elinewidth=1, color=(0,0,0,0.3), zorder=-1)
        
        ax.legend(fontsize=7, loc='lower right', framealpha=0.2, borderpad=0.3)
        
        upe = rpar[t]
        upeph = rpar[t + '_pH']
        if t == 'white':
            parlab = ("A: {:.3f}\n".format(upeph[0]) +
                    "B: {:.3f}\n".format(upeph[1]) +
                    "$C_1$" + ": {:.1f}\n".format(upeph[2]) +
        #               "$C_2 (CO3)$" + ": {:.3f}\n".format(upeph[3]) + 
                    "$C_2$" + ": {:.3f}$\\times 10^6$\n".format(upeph[3] * 1e-6) + 
                    "D: {:.3f}\n".format(upeph[4]))
        else:
            parlab = ("P: {:.3f}\n".format(upeph[0]) +
                    "$C_1$" + ": {:.1f}\n".format(upeph[1]) +
        #               "$C_2 (CO3)$" + ": {:.3f}\n".format(upeph[3]) + 
                    "$C_2$" + ": {:.3f}$\\times 10^6$\n".format(upeph[2] * 1e-6) + 
                    "D: {:.3f}\n".format(upeph[3]))
        y = 0.76
        # if i == 1:
        #     y = 0.76
        # else:
        #     y = 0.84
        ax.text(.03, y, parlab, va='top', ha='left', fontsize=6, color=(.3, .3, .3), transform=ax.transAxes)
        
    # Gray data
    ax = axs[2]
    ax.scatter(gray.loc[:, 'G. ruber w Mg/Ca [mmol/mol]'], nom(rpred['gray_pH']), 
            marker='.', label='Gray et al. (2018)', s=15, **tstyle['white'], alpha=0.6)
    ax.legend(fontsize=7, loc='lower right', framealpha=0.2, borderpad=0.3)

    upe = rpar['gray']
    upeph = rpar['gray_pH']
    parlab = ("A: {:.3f}\n".format(upeph[0]) +
            "B: {:.3f}\n".format(upeph[1]) +
            "$C_1$" + ": {:.1f}\n".format(upeph[2]) +
    #           "$C_2 (CO3)$" + ": {:.3f}\n".format(upeph[3]) + 
            "$C_2$" + ": {:.3f}$\\times 10^6$\n".format(upeph[3] * 1e-6) + 
            "D: {:.3f}\n".format(upeph[4]))
    ax.text(.05, .84, parlab, va='top', ha='left', fontsize=6, color=(.3, .3, .3), transform=ax.transAxes)

    desc = ['culture - white', 'culture - white', 'sed. trap - white', 'culture - pink']
    # axis labels
    for i, ax in enumerate(axs):
        ax.set_xlim(lim)
        ax.set_ylim(lim[0], lim[-1] * 1.15)
        ax.plot(lim, lim, zorder=-1, ls='dashed', color=(0,0,0,0.4))
        
        if ax.is_last_row():
            ax.set_xlabel('Measured $Mg/Ca_{\it{G. ruber}}$ (mmol $mol^{-1}$)')
        if ax.is_first_col():
            ax.set_ylabel('Predicted $Mg/Ca_{\it{G. ruber}}$ (mmol $mol^{-1}$)')
        ax.text(.05, .95, labels[i], ha='left', va='top', transform=ax.transAxes, weight='bold',
                fontsize=14, color=(.3, .3, .3), zorder=2)
        ax.text(.125, .935, desc[i], ha='left', va='top', transform=ax.transAxes, weight='bold', style='italic',
                fontsize=8, color=(.6, .6, .6), zorder=2)

    ax = axs[1]
    eqnlab = '$Mg/Ca_{foram} = Mg/Ca_{SW}^A\ B\ e^{(C_1\ [Ca]_{SW}\ + C_2\ [CO_3]_{SW} + D) T}$\n'
    ax.text(.03, .86, eqnlab, va='top', ha='left', fontsize=7, color=(.15, .15, .15), transform=ax.transAxes)
    
    ax = axs[3]
    eqnlab = '$Mg/Ca_{foram} = P\ e^{(C_1\ [Ca]_{SW}\ + C_2\ [CO_3]_{SW} + D) T}$\n'
    ax.text(.03, .86, eqnlab, va='top', ha='left', fontsize=7, color=(.15, .15, .15), transform=ax.transAxes)

        
    fig.tight_layout(h_pad=0.01)

    for i, p in zip([1,3], ['white', 'pink']):
        upe = rpar[p]
        
    # CO3 colour scale
    pos = axs[0].get_position()
    cbax = fig.add_axes([pos.x1 - pos.width * 0.45, pos.y0 + pos.height * 0.15,
                        pos.width * 0.4, pos.height * 0.05])
    fig.colorbar(cb, cax=cbax, orientation='horizontal', ticks=[7.5, 8, 8.5], label='$pH_{Total}$')
    cbax.xaxis.set_label_position('top')

    # ticks
    for ax in axs:
        ax.xaxis.set_major_locator(ticker.MultipleLocator(2))
        ax.yaxis.set_major_locator(ticker.MultipleLocator(2))
        ax.xaxis.set_minor_locator(ticker.MultipleLocator(.5))
        ax.yaxis.set_minor_locator(ticker.MultipleLocator(.5))
    

    return fig, axs

# code for figure 6 is in online supplement

def fig7(xn, yn_ca, yn_mgca, ref_ca, ref_mg, yr_mgca, xr_mgca, MgCa_CC, mx, my, dca):
    m = ['o', 'v', '^', '<', '>', '_', '|', '+', 'x', 's', 'D', '*']
    rs = np.unique(np.concatenate([ref_mg, ref_ca]))
    mdict = {r: m for r, m in zip(rs, m[:len(rs)])}

    w = 88 / 25.4
    size = (1.8 * w, w * 2.2)

    cm = plt.cm.Blues_r
    cm.set_bad(  (.8, .8, .8))
    cm.set_under((.8, .8, .8))
    cm.set_over( (.8, .8, .8))

    cmg = plt.cm.Reds
    cmg.set_bad(  (.8, .8, .8))
    cmg.set_under((.8, .8, .8))
    cmg.set_over( (.8, .8, .8))


    fig = plt.figure(figsize=size)

    cax = fig.add_axes([.14,.78,.8 / 1.8,.19])
    mgcax = fig.add_axes([.14,.58,.8 / 1.8,.19])
    tsx = fig.add_axes([.14,.38,.8 / 1.8,.19])

    cbx = fig.add_axes([.14 + .8 / 1.8 + 0.01,.385,.015,.18])

    # labels
    cax.text(.02, .95, 'A', ha='left', va='top', transform=cax.transAxes, weight='bold',
            fontsize=14, color=(.3, .3, .3), zorder=2)
    for ax, label in zip([mgcax, tsx], ['B', 'C']):
        ax.text(.98, .95, label, ha='right', va='top', transform=ax.transAxes, weight='bold',
                fontsize=14, color=(.3, .3, .3), zorder=2)

    # [Ca]
    cax.plot(xn, nom(yn_ca), c='k', alpha=0.8)
    cax.fill_between(xn, nom(yn_ca) - err(yn_ca) * 1.96, nom(yn_ca) + err(yn_ca) * 1.96, 
                    zorder=-1, alpha=0.2, color='k')

    ca_legend = []
    for i, r in dca.iterrows():
        R,X,L,T,Y,B = r.loc[idx[('Age (Myr)', 'Casw (mM)'),:]].values
        
        if all(~r.loc[idx[('Age (Myr)', 'Casw (mM)'),:]].isnull()):
            p = mpl.patches.Polygon([[L,B],[L,T],[R,T],[R,B]], alpha=0.2, color='k')
            cax.add_patch(p)
        ref = r.loc[idx['Reference','Reference']]
        if (np.nanmin(Y) < 120) and ref not in ca_legend:
            label = ref
            ca_legend.append(ref)
        else:
            label = "_"
        cax.scatter(X, Y, marker=mdict[ref],
                    c='k', zorder=0, alpha=0.2, s=15, label=label)
            
    cax.legend(loc='upper left', bbox_to_anchor=(1.01,1), fontsize=7)

    cax.set_ylabel('$[Ca]_{SW}$ (mmol $kg^{-1}$)')
        
    # Mg/Ca
    mgcax.plot(xn, yn_mgca['mean'], c='k', alpha=0.8)
    mgcax.fill_between(xn, yn_mgca['mean'] - yn_mgca['std'], 
                    yn_mgca['mean'] + yn_mgca['std'], zorder=-1, alpha=0.2, color='k')

    mg_legend = []
    for r in np.unique(ref_mg):
        ind = ref_mg == r
        if (np.nanmin(yr_mgca[ind]) < 120) and r not in mg_legend:
            label = r
        else:
            label = "_"

        mgcax.scatter(xr_mgca[ind], yr_mgca[ind], c='k', zorder=0, alpha=0.2, s=15, label=label, marker=mdict[r])
    mgcax.legend(loc='upper left', bbox_to_anchor=(1.01,1), fontsize=7)

    mgcax.set_ylabel('$Mg/Ca_{SW}$ (mol $mol^{-1}$)')

    # temperature sensitivity axis
    vmin = np.nanmin(MgCa_CC)
    vmax = np.nanmax(MgCa_CC)

    cmg = tsx.pcolormesh(mx, my, MgCa_CC, cmap=cm, zorder=-2, vmin=vmin, vmax=vmax)
    contours = [2, 4, 6, 8, 10, 12]
    csl = tsx.contour(mx, my, MgCa_CC, contours, colors=[(0,0,0,0.3)], linewidths=1)
    contour_labels(csl, tsx, fmt='%.0f', hshift=-.8)
    fig.colorbar(cmg, cbx, orientation='vertical', label='Mg/Ca (mmol $mol^{-1}$)', ticks=contours)

    tsx.set_ylabel('Temperature ($^{\circ}C)$')
    tsx.set_xlabel('Time (Mya)')

    for ax in [cax, mgcax, tsx]:
        ax.set_xlim(0,120)

    for ax in [cax, mgcax]:
        ax.set_xticklabels([])

    return fig

def fig8(params):
    # Conditions
    Mg = 50e-3
    Ca = 10.2e-3
    MgCa = Mg/Ca
    DIC = 2000e-6
    MgCa_CC = 5

    Ca_Paleocene = 20e-3
    MgCa_Paleocene = 1.5
    MgCa_CC_Paleocene = 5
    DIC_Paleocene = 2000e-6

    c_Paleocene = 'C1'
    c_Modern = 'C0'

    fig, axs = plt.subplots(2, 2, figsize=[6, 5])

    axs = axs.flatten()

    labels = 'ABCDEFG'

    for mgcaf, ls in zip((5, 7), ('solid', 'dashed')):
        MgCa_CC = mgcaf
        MgCa_CC_Paleocene = mgcaf
        
        # Mg/Casw uncertainty
        ax = axs[0]
        uMgCa = un.unumpy.uarray(MgCa, np.linspace(0, .5, 50))
        Tpred = T_fn((uMgCa, Ca, DIC, MgCa_CC), *params)
        ax.plot(err(uMgCa),
                err(Tpred), c=c_Modern, ls=ls)

        uMgCa_Paleocene = un.unumpy.uarray(MgCa_Paleocene, np.linspace(0, .5, 50))
        Tpred_Paleocene = T_fn((uMgCa_Paleocene, Ca_Paleocene, DIC_Paleocene, MgCa_CC_Paleocene), *params)

        ax.plot(err(uMgCa_Paleocene),
                err(Tpred_Paleocene), c=c_Paleocene, ls=ls)
        
        # Ca uncertainty
        ax = axs[1]
        uCa = un.unumpy.uarray(Ca, np.linspace(0, 2e-3, 50))
        Tpred = T_fn((MgCa, uCa, DIC, MgCa_CC), *params)

        ax.plot(err(uCa) * 1e3,
                err(Tpred), c=c_Modern, ls=ls)
        

        uCa_Paleocene = un.unumpy.uarray(Ca_Paleocene, np.linspace(0, 2e-3, 50))
        Tpred_Paleocene = T_fn((MgCa_Paleocene, uCa_Paleocene, DIC_Paleocene, MgCa_CC_Paleocene), *params)
        ax.plot(err(uCa_Paleocene) * 1e3,
                err(Tpred_Paleocene), c=c_Paleocene, ls=ls)

        # DIC uncertainty
        ax = axs[2]
        uDIC = un.unumpy.uarray(DIC, np.linspace(0, 600e-6, 50))
        Tpred = T_fn((MgCa, Ca, uDIC, MgCa_CC), *params)

        ax.plot(err(uDIC) * 1e6,
                err(Tpred), c=c_Modern, ls=ls)

        uDIC_Paleocene = un.unumpy.uarray(DIC_Paleocene, np.linspace(0, 600e-6, 50))
        Tpred_Paleocene = T_fn((MgCa_Paleocene, Ca_Paleocene, uDIC_Paleocene, MgCa_CC_Paleocene), *params)
        ax.plot(err(uDIC_Paleocene) * 1e6,
                err(Tpred_Paleocene), c=c_Paleocene, ls=ls)

        # Mg/Ca foram uncertainty
        ax = axs[3]
        uMgCa_CC = un.unumpy.uarray(MgCa_CC, np.linspace(0, 0.4, 50))
        Tpred = T_fn((MgCa, Ca, DIC, uMgCa_CC), *params)

        ax.plot(err(uMgCa_CC),
                err(Tpred), c=c_Modern, ls=ls)

        uMgCa_CC_Paleocene = un.unumpy.uarray(MgCa_CC_Paleocene, np.linspace(0, .4, 50))
        Tpred_Paleocene = T_fn((MgCa_Paleocene, Ca_Paleocene, DIC_Paleocene, uMgCa_CC_Paleocene), *params)
        ax.plot(err(uMgCa_CC),
                err(Tpred_Paleocene), c=c_Paleocene, ls=ls)

    axs[0].set_xlabel('$\sigma\ Mg/Ca_{SW}\ (mol\ mol^{-1})$')
    axs[1].set_xlabel('$\sigma\ [Ca]_{SW}\ (mmol\ kg^{-1})$')
    axs[2].set_xlabel('$\sigma\ [DIC]_{SW}\ (\mu mol\ kg^{-1})$')
    axs[3].set_xlabel('$\sigma\ Mg/Ca_{\it{O. universa}}\ (mmol\ mol^{-1})$')

    #legend
    axs[1].plot([], [], c=c_Modern, label='Modern')
    axs[1].plot([], [], c=c_Paleocene, label='Paleocene')
    axs[1].plot([], [], c=(0,0,0,0.7), label='$Mg/Ca_{\it{O. universa}}$ = 5 mmol $mol^{-1}$', ls='solid')
    axs[1].plot([], [], c=(0,0,0,0.7), label='$Mg/Ca_{\it{O. universa}}$ = 7 mmol $mol^{-1}$', ls='dashed')
    axs[1].legend(fontsize=8, framealpha=0.6)

    axs[0].set_ylabel('$\sigma\ T\ (^{\circ}C)$')
    axs[1].set_ylim(axs[0].get_ylim())
    axs[1].set_yticklabels([])
    axs[2].set_ylabel('$\sigma\ T\ (^{\circ}C)$')
    axs[2].set_ylim(axs[3].get_ylim())
    axs[3].set_yticklabels([])

    for i, ax in enumerate(axs):
        ax.text(.03, .96, labels[i], ha='left', va='top', transform=ax.transAxes, weight='bold',
                fontsize=14, color=(.3, .3, .3), zorder=2)

    fig.tight_layout()

    return fig, axs
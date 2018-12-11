import pandas as pd
import numpy as np
import cbsyst as cb
from tqdm import tqdm
from mg_funks import phreeq as ph

def mg_data(excel_file, species=None):
    """
    Loads Mg/Ca data.
    
    Parameters
    ----------
    excel_file : str
        Path to excel file containing data.
    species : str
        A string containing either "Genus species", "Genus" or "species",
        specifying which species to load.
    """
    d = pd.read_excel(excel_file)
    ind = np.ones(d.shape[0], dtype=bool)
    
    if species is not None:
        sspl = species.split(' ')
        if len(sspl) == 2:
            genus = sspl[0]
            speci = sspl[1]
            ind = ind & (d.loc[:, 'Genus'] == genus) & (d.loc[:, 'Species'] == speci)
        elif len(sspl) == 1:
            genus = sspl[0]
            speci = sspl[0]
            ind = ind & ((d.loc[:, 'Genus'] == genus) | (d.loc[:, 'Species'] == speci))
        else:
            vspec = '\n    - ' + '\n    - '.join((d.loc[:, 'Genus'] + ' ' + d.loc[:, 'Species']).unique())
            raise ValueError('Dont understand species="{}". Please use a valid species name:'.format(species) + vspec)

    d.loc[:, 'who'] = [s[0] for s in d.Reference.str.split(',')]

    return d.loc[ind, :].apply(pd.to_numeric, errors='ignore')

def calc_myami_Csys(dat):
    cols = ['pHtot_out', 'DIC_out', 'HCO3_out', 'CO3_out', 'TA_out']

    dic = dat.loc[:, 'DIC'].copy()
    # pHNBS and DIC
    ind = ~np.isnan(dat.loc[:, 'pHNBS']) & ~np.isnan(dic)
    if any(ind):
        cs = cb.Csys(pHNBS=dat.loc[ind, 'pHNBS'], DIC=dic.loc[ind], BT=dat.loc[ind, 'B umol/kg'],
                     Ca=dat.loc[ind, '[Ca]sw'] * 1e-3, Mg=dat.loc[ind, '[Mg]sw'] * 1e-3, T_in=25., T_out=dat.loc[ind, 'Temp'],
                     S_in=dat.loc[ind, 'Salinity'],)

        for c in cols:
            dat.loc[ind, c.replace('_out', '')] = cs[c]
            dat.loc[ind, 'omega'] = cs.CO3_out * 1e-6 * cs.Ca / cs.Ks.KspC
    
    # pHTOT and DIC
    ind = ~np.isnan(dat.loc[:, 'pHTOTAL']) & ~np.isnan(dic)
    if any(ind):
        cs = cb.Csys(pHtot=dat.loc[ind, 'pHTOTAL'], DIC=dic.loc[ind], BT=dat.loc[ind, 'B umol/kg'],
                     Ca=dat.loc[ind, '[Ca]sw'] * 1e-3, Mg=dat.loc[ind, '[Mg]sw'] * 1e-3, T_in=25., T_out=dat.loc[ind, 'Temp'],
                     S_in=dat.loc[ind, 'Salinity'],)

        for c in cols:
            dat.loc[ind, c.replace('_out', '')] = cs[c]
            dat.loc[ind, 'omega'] = cs.CO3_out * 1e-6 * cs.Ca / cs.Ks.KspC

    # pHTOT and TA
    ind = ~np.isnan(dat.loc[:, 'pHTOTAL']) & ~np.isnan(dat.loc[:, 'Alk'])
    if any(ind):
        cs = cb.Csys(pHtot=dat.loc[ind, 'pHTOTAL'], TA=dat.loc[ind, 'Alk'], BT=dat.loc[ind, 'B umol/kg'],
                     Ca=dat.loc[ind, '[Ca]sw'] * 1e-3, Mg=dat.loc[ind, '[Mg]sw'] * 1e-3, T_in=25., T_out=dat.loc[ind, 'Temp'],
                     S_in=dat.loc[ind, 'Salinity'],)

        for c in cols:
            dat.loc[ind, c.replace('_out', '')] = cs[c]
            dat.loc[ind, 'omega'] = cs.CO3_out * 1e-6 * cs.Ca / cs.Ks.KspC

    # pHNBS and TA
    ind = ~np.isnan(dat.loc[:, 'pHNBS']) & ~np.isnan(dat.loc[:, 'Alk'])
    if any(ind):
        cs = cb.Csys(pHNBS=dat.loc[ind, 'pHNBS'], TA=dat.loc[ind, 'Alk'], BT=dat.loc[ind, 'B umol/kg'],
                     Ca=dat.loc[ind, '[Ca]sw'] * 1e-3, Mg=dat.loc[ind, '[Mg]sw'] * 1e-3, T_in=25., T_out=dat.loc[ind, 'Temp'],
                     S_in=dat.loc[ind, 'Salinity'],)

        for c in cols:
            dat.loc[ind, c.replace('_out', '')] = cs[c]
            dat.loc[ind, 'omega'] = cs.CO3_out * 1e-6 * cs.Ca / cs.Ks.KspC

    return dat

def calc_pitzer_Csys(dat, 
                     database='/home/oscar/phreeqc/iphreeqc-3.3.9-11951/database/pitzer.dat', 
                     phreeq_path='/usr/local/lib/libiphreeqc.so'):
    """
    Function to calculate C specitation, corrected for [Mg] and [Ca] using PHREEQC.

    Parameters
    ----------
    d : pandas.DataFrame
        Output of mg_data load function.
    database : str
        Path to PITZER database file.
    phreeq_path : str
        Path to PHREEQC executable.

    """
    # remove data without [Mg] and [Ca]
    dat = dat.loc[~dat.loc[:,['[Mg]sw', '[Ca]sw']].isnull().any(1)]
    # create working copy
    d = dat.copy()
    # rename columns
    d.columns = pd.MultiIndex.from_product([['Measured'], d.columns])
    
    # define standard seawater composition (mol/kg SW)
    # ================================================
    sw = {'pH': 8.2,
          'Temp': 25,
          'units': 'mol/kgw',
          'density': 1.026,
          'Ca': 10.28e-3,
          'Cl': 0.54586,
          'K': 10.21e-3,
          'Mg': 52.82e-3,
          'Na': 0.46906,
          'S(6)': 0.02824,
          'B': 416e-6,
          'C': 2.3e-3,
          'Br': 0.84e-3,
          'F': 0.7e-4}
    
    # calculate ambient Ks at 25C (measurement temperature)
    # =====================================================
    ks = {}

    TempC = 25
    TempK = TempC + 273.15
    Sal = d.loc[:, ('Measured', 'Salinity')]

    par = cb.MyAMI_V2.start_params
    fns = cb.MyAMI_V2.fn_dict

    for k in ['K0', 'K1', 'K2', 'KB', 'KW', 'KspC', 'KspA', 'KSO4']:
        ks[k] = fns[k]((TempK, Sal), *par[k])

    # calculate other Ks
    ks.update(cb.non_MyAMI_constants.calc_KF(TempC, Sal))
    ks.update(cb.non_MyAMI_constants.calc_KPs(TempC, Sal))
    ks.update(cb.non_MyAMI_constants.calc_KSi(TempC, Sal))
    
    # calculate correction factors at 25C
    # ===================================
    fKs = {}
    isw = sw.copy()

    salsc = ['Cl', 'K', 'Na', 'S(6)', 'Br', 'F']

    for i, v in tqdm(d.iterrows(), total=d.shape[0],
                     desc='Calculating 25C Corrections'):
        v = v.Measured
        isw.update({'Mg': v.loc['[Mg]sw'] * 1e-3,
                    'Ca': v.loc['[Ca]sw'] * 1e-3,
                    'Temp': 25,
                    'B': v.loc['B umol/kg'] * 1e-6})

        sal_f = v.Salinity / 35.
        isw.update({k: sw[k] * sal_f for k in salsc})

        fks = ph.calc_pitzer_fKs(sw, isw, database=database, phreeq_path=phreeq_path)

        for k, f in fks.items():
            if k not in fKs:
                fKs[k] = []
            fKs[k].append(f)
    
    # apply correction factors to 25C Ks
    # ==================================
    cKs = ks.copy()
    for k, f in fKs.items():
        cKs[k] = ks[k] * f
    
    # convert all pH scales to Total
    # ==============================
    pHs = cb.helpers.calc_pH_scales(pHtot=None, pHfree=None, pHsws=None,
                                    pHNBS=d.loc[:, ('Measured', 'pHNBS')],
                                    TS=cb.calc_TS(Sal), TF=cb.calc_TF(Sal),
                                    TempK=TempK, Sal=Sal, Ks=cb.Bunch(cKs))

    ind = d.loc[:, ('Measured', 'pHTOTAL')].isnull()
    d.loc[ind, ('Measured', 'pHTOTAL')] = pHs['pHtot'][ind]
    
    # Calculate DIC/Alk at MEASURED conditions (25C)
    # ==============================================
    # calculate carbon system using 25C Ks, pH and either DIC or Alk
    carb = cb.Csys(d.loc[:, ('Measured', 'pHTOTAL')], TA=d.loc[:, ('Measured', 'Alk')], 
                   BT=d.loc[:, ('Measured', 'B umol/kg')], Ks=cKs)
    d.loc[:, ('pitzer_25C', 'DIC')] = carb.DIC

    carb = cb.Csys(d.loc[:, ('Measured', 'pHTOTAL')], DIC=d.loc[:, ('Measured', 'DIC')], 
                   BT=d.loc[:, ('Measured', 'B umol/kg')], Ks=cKs)
    d.loc[:, ('pitzer_25C', 'Alk')] = carb.TA

    # whichever of DIC or Alk wasn't calculated, transfer the
    # measured value to the pitzer_25C column for in-situ condition calculation
    dicnull = d.loc[:, ('pitzer_25C', 'DIC')].isnull()
    d.loc[dicnull, ('pitzer_25C', 'DIC')] = d.loc[dicnull, ('Measured', 'DIC')]

    alknull = d.loc[:, ('pitzer_25C', 'Alk')].isnull()
    d.loc[alknull, ('pitzer_25C', 'Alk')] = d.loc[alknull, ('Measured', 'Alk')]
    
    
    # calculate empirical Ks for ambient Mg and Ca, and experiment Sal and Temp
    # =========================================================================
    eKs = cb.calc_Ks(d.loc[:, ('Measured', 'Temp')], d.loc[:, ('Measured', 'Salinity')], 0, sw['Mg'], sw['Ca'],
                     cb.calc_TS(d.loc[:, ('Measured', 'Salinity')]), cb.calc_TF(d.loc[:, ('Measured', 'Salinity')]))
    
    # calculate correction factors at experimental temperature
    # ========================================================
    fKs = {}
    isw = sw.copy()

    salsc = ['Cl', 'K', 'Na', 'S(6)', 'Br', 'F']

    for i, v in tqdm(d.iterrows(), total=d.shape[0],
                     desc='Calculating Experimental T Corrections'):
        v = v.Measured
        isw.update({'Mg': v.loc['[Mg]sw'] * 1e-3,
                    'Ca': v.loc['[Ca]sw'] * 1e-3,
                    'Temp': v.Temp,
                    'B': v.loc['B umol/kg'] * 1e-6})

        sal_f = v.Salinity / 35.
        isw.update({k: sw[k] * sal_f for k in salsc})

        fks = ph.calc_pitzer_fKs(sw, isw, database=database, phreeq_path=phreeq_path)

        for k, f in fks.items():
            if k not in fKs:
                fKs[k] = []
            fKs[k].append(f)
    
    
    # apply correction factors
    for k, f in fKs.items():
        eKs[k] = ks[k] * f

    # calculate carbon system at experimental conditions using pitzer-calculated DIC/Alk pairs
    carb = cb.Csys(DIC=d.loc[:, ('pitzer_25C', 'DIC')], TA=d.loc[:, ('pitzer_25C', 'Alk')],
                   BT=d.loc[:, ('Measured', 'B umol/kg')], Ks=cKs)

    # store outputs!
    for c in ['DIC', 'TA', 'CO3', 'HCO3', 'CO2', 'pHtot']:
        d.loc[:, ('pitzer_expT', c)] = carb[c]
        dat.loc[:, c] = carb[c]

    return dat


def orbulina(data_file):
    """
    Load all Orbulina universa data.
    """

    # Kate (Ca and Mg experiments and T)
    kd = pd.read_excel(data_file, sheet_name='Kate', skiprows=1)
    kd = kd.iloc[1:,:].dropna()
    kd.loc[:,'who'] = 'Kate'
    
    # Kate 2 (carbonate experiments)
    kd2 = pd.read_excel(data_file, sheet_name='Kate2', skiprows=1)
    kd2 = kd2.iloc[1:,].dropna()
    kd2.loc[:,'who'] = 'Kate2'
    # exclude low CO3 points
    # kd2 = kd2.loc[kd2.loc[:,'CO32-'] > 50,:]

    # Steve Ca and Temp and Mg
    sd = pd.read_excel(data_file, sheet_name='Steve', skiprows=1)
    sd = sd.iloc[1:, :]
    sd.loc[:,'Mg/Caf'] /= 1000
    sd.loc[:,'Mg/Caf 2se'] /= 1000
    sd.loc[:,'who'] = 'Steve'

    # Ann
    ar = pd.read_excel(data_file, sheet_name='Russell', skiprows=1)
    ar = ar.iloc[1:]
    ar.loc[:,'who'] = 'Ann'
    ar.loc[:,'Mg/Caf'] /= 1000
    ar.loc[:,'Mg/Caf 2se'] /= 1000

    # exclude low CO3 point BECAUSE IT'S WEIRD!!!
    # ar = ar.loc[ar.loc[:,'CO32-'] > 100, :]

    # Kat Allen
    ka = pd.read_excel(data_file, sheet_name='Allen', skiprows=1)
    ka = ka.iloc[1:]
    ka.loc[:,'who'] = 'Allen'
    ka.loc[:,'Mg/Caf'] /= 1000
    ka.loc[:,'Mg/Caf 2se'] /= 1000

    # Can't use Kat's data - no coherent pattern with T!

    # Barbel
    hd = pd.read_excel(data_file, sheet_name='Honisch', skiprows=1)
    hd = hd.iloc[1:,:].dropna()
    hd.loc[:,'who'] = 'Barbel'

    # Lozza
    zd = pd.read_excel(data_file, sheet_name='Loz', skiprows=1)
    zd = zd.iloc[2:,:]
    zd.loc[:,'who'] = 'Lozza'

    # combine data
    dat = pd.concat([kd2, ar, kd, sd, hd, zd, ka])

    # data exclusion
    ind = dat.pH > 7.5
    print("Note: Excluded Data:")
    print(dat.loc[~ind])

    dat = dat.loc[ind]

    for c in ['[Ca]sw', 'Mg/Caf', 'Mg/Casw', '[Mg]sw']:
        dat.loc[:, c] = dat.loc[:, c].astype(float)
    
    # placeholder: replace missing errors with median uncertainty
    dat.loc[:, 'Mg/Caf 2se'].fillna(dat.loc[:, 'Mg/Caf 2se'].quantile(0.5), inplace=True)

    return dat  
#!/usr/bin/env python

"""This script defines functions to equilibrium simulation of gasification 
processes. It uses some predefined functions from Cantera package.

@author = Rodolfo Rodrigues
@contact = rodolfo.rodrigues@ufsm.br
@data = April, 2012, rev.: June, 2013 (adapted to use cython Cantera)

TODO: 2022/05/02 Create function to convert between air, o2, ER and between steam, SR.
"""

#==============================================================================
# import libraries/files
#==============================================================================
import pp as ppold
import pp2 as pp
import feedstock as fsold
import feedstock2 as fs
import fuel as fu
import cantera as ct
import numpy as np
import scipy.optimize as opt
import csv

#==============================================================================
# predefine parameters
#==============================================================================
R = ct.gas_constant  # 8314.4621 Pa*m^3/K/kmol
Tn = 273.15  # K
Pn = ct.one_atm  # 101315 Pa
zero = np.zeros(1)
one = np.ones(1)

#==============================================================================
# special functions
#==============================================================================

def getFeed(fuelMix, moist=0.0, air=0.0, steam=0.0):
    '''
    This function creates a mixture of phases to denote the fuel.
    The fuel is composed as a mixture of char, gas, ash, and moisture phases.

    Parameters
    ----------
    fuelMix : Cantera 'Mixture' object
        Object containing the mole amount of each species in the dry fuel.
    moist : float
        Mass fraction of moisture fuel [kg/kg] (default value is zero)
    air : float
        Mass amount of air [kg] (default value is zero)
    steam : float
        Mass amount of steam [kg] (default value is zero)

    Returns
    -------
    feed : Cantera 'Mixture' object
        Object representing the mixture of phases in the feedstock.
    '''

    feed = pp.mix()

    mw = np.fromiter(pp.Mw.values(), dtype=float) # molecular weights
    fuelMoles = fuelMix.species_moles # moles for each species
    fuelMass = fuelMoles * mw / 1000 # fuel mass in kg for each species
    totalFuelMass = np.sum(fuelMass) # total fuel mass in kg

    moistMass = moist * totalFuelMass # mass of d.b. moisture in kg
    moistMoles = moistMass * 1000 / pp.Mw['H2O'] # moles of d.b. moisture
    
    steamMoles = steam * 1000 / pp.Mw['H2O'] # moles of steam

    H2OMoles = moistMoles + steamMoles # total moles of water

    airO2Moles = 1000*0.23211606*air/pp.Mw['O2']
    airN2Moles = 1000*0.75507754*air/pp.Mw['N2']
    airArMoles = 1000*0.01280640*air/pp.Mw['Ar']

    #pureO2Moles = O2 / pp.Mw['O2'] # moles of pure O2

    O2Moles = airO2Moles #+ pureO2Moles total moles of O2

    feedMoles = np.zeros(len(feed.species_names))
    feedMoles += fuelMoles
    feedMoles[pp.i['H2O']] += H2OMoles
    feedMoles[pp.i['O2']] += O2Moles
    feedMoles[pp.i['N2']] += airN2Moles
    feedMoles[pp.i['Ar']] += airArMoles
    feed.species_moles = feedMoles

    return feed

def isotGasification(fuelID, fuelMass, moisture, T=1273.15, P=ct.one_atm, 
                    oxi = 0.5, steam=0.0, oType='ER', sType='steam',
                    species=['C(gr)','N2','O2','H2','CO','CH4','CO2','H2O']):
    '''
    Isothermal gasification calculation for a single fuel in a given condition.

    Parameters
    ----------
    fuelID : str
        ID of fuel as given by the database (fuels.csv)
    fuelMass : float
        The fuel mass [kg]
    moisture : float
        The moisture mass fraction, dry base [kg/kg]
    T : float
        Temperature [K]
    P : float
        Pressure [atm]
    oxi : float
        Oxidizer value [kg air] [kg O2] [kg/kg]
    steam : float
        Steam value [kg] [kg/kg]
    oType : str
        Oxidizer type, 'air', 'O2' or 'ER'
    sType : str
        Steam type, 'steam' or 'SR'
    species : list
        List of species to include in the calculation.


    Returns
    -------
    outlet : Cantera 'Mixture' object
        Object representing the mixture at equilibrium.
    '''
    # Create fuel mix
    fuelMix = fs.getFuelMix(fuelID, fuelMass)

    # Verify type of oxidizer given
    if oType == 'air':
        air = oxi
        O2 = 0
        ER = fs.airtoER(fuelMix, oxi)
    elif oType == 'O2':
        air = 0
        O2 = oxi
        ER = fs.airtoER(fuelMix, 0, oxi) # TODO: Can't use pure O2
    elif oType == 'ER':
        air = fs.ERtoair(fuelMix, oxi)
        O2 = 0
        ER = oxi
    else:
        raise ValueError('Invalid oxidizer type')
    
    # Verify type of steam given
    if sType == 'steam':
        stm = steam
        SR = fs.SRtosteam(fuelMix, steam)
    elif sType == 'SR':
        stm = fs.steamtoSR(fuelMix, steam)
        SR = steam
    else:
        raise ValueError('Invalid steam type')

    # Create feed
    feed = getFeed(fuelMix, moisture, air, steam)

    # Calculate equilibrium
    feed.T = T
    feed.P = P
    feed.equilibrate('TP')

    return feed

def isotCogasification(fuel1, fuel2, fuel1Mass, blend, moisture, T=1273.15, 
                    P=ct.one_atm, oxi = 0.5, steam=0.0, oType='ER', sType='steam',
                    species=['C(gr)','N2','O2','H2','CO','CH4','CO2','H2O']):
    '''
    Isothermal gasification calculation for a blend of 2 fuels in a given condition.

    Parameters
    ----------
    fuel1 : str
        ID of fuel #1 as given by the database (fuels.csv)
    fuel2 : str
        ID of fuel #2 as given by the database (fuels.csv)
    fuel1Mass : float
        Mass of fuel #1 [kg]
    blend : float
        Fuel #2 to total fuel mass ratio [kg]
    moisture : float
        The moisture mass fraction, dry base [kg/kg]
    T : float
        Temperature [K]
    P : float
        Pressure [atm]
    oxi : float
        Oxidizer value [kg air] [kg O2] [kg/kg]
    steam : float
        Steam value [kg] [kg/kg]
    oType : str
        Oxidizer type, 'air', 'O2' or 'ER'
    sType : str
        Steam type, 'steam' or 'SR'
    species : list
        List of species to include in the calculation.


    Returns
    -------
    outlet : Cantera 'Mixture' object
        Object representing the mixture at equilibrium.
    '''
    # Create each fuel mixture
    fuel1Mix = fs.getFuelMix(fuel1, fuel1Mass)
    fuel2Mass = fuel1Mass * blend # WRONG
    fuel2Mix = fs.getFuelMix(fuel2, fuel2Mass)

    # Create fuel blend
    fuelBlend = fs.blend(fuel1Mix, fuel2Mix)

    # Verify type of oxidizer given
    if oType == 'air':
        air = oxi
        O2 = 0
        ER = fs.airtoER(fuelBlend, oxi)
    elif oType == 'O2':
        air = 0
        O2 = oxi
        ER = fs.airtoER(fuelBlend, 0, oxi) # TODO: Can't use pure O2
    elif oType == 'ER':
        air = fs.ERtoair(fuelBlend, oxi)
        O2 = 0
        ER = oxi
    else:
        raise ValueError('Invalid oxidizer type')
    
        # Verify type of steam given
    if sType == 'steam':
        stm = steam
        SR = fs.SRtosteam(fuelBlend, steam)
    elif sType == 'SR':
        stm = fs.steamtoSR(fuelBlend, steam)
        SR = steam
    else:
        raise ValueError('Invalid steam type')

    # Create feed
    feed = getFeed(fuelBlend, moisture, air, steam)

    # Calculate equilibrium
    feed.T = T
    feed.P = P
    feed.equilibrate('TP')

    return feed


            
# def coprocessing(self, fuel_id, blend, moisture, T, P=1.0,
#                  air=0, O2=0, ER=0.4, steam=0, SR=0,
#                  small=None, db='n', normalized='n', format_='%',
#                  species=['C(gr)','N2','O2','H2','CO','CH4','CO2','H2O']):
#     """
#     Cogasification calculations for binary blends of fuels.
    
#     Parameters
#     ----------
#     self : ndarray
#         Mass fraction of fuel #1 [kg/kg]
#     fuel_id : list of strings
#         List of ID fuel
#     blend : float|ndarray
#         Fuel #1 to fuel #2 ratio [kg/kg]
#     moisture : float|ndarray
#         Moisture mass fraction [kg/kg]
#     T : float|ndarray
#         Temperature [degC]
#     P : float|ndarray
#         Pressure [atm] (default is 1.0)
#     air : float|ndarray
#         Air amount [kg] (default is zero)
#     O2 : float|ndarray
#         O2 amount [kg] (default is zero)
#     ER : float|ndarray
#         Equivalence ratio [kmol/kmol]
#     steam : float|ndarray
#         Steam amount [kg] (default is zero)
#     SR : float|ndarray
#         Steam to carbon ratio [kmol/kmol] (default is zero)
#         basis: 1 kg coal-biomass blend, d.b.
#     small : float
#         Smallest number to report as a fraction value (default is None)
#     db : string
#         Dry basis composition ('y') or web basis composition ('n') (default
#         is 'n')
#     normalized : string
#         Normalized compostion ('y') or overall composition ('n') (default
#         is 'n')
#     format_ : string
#         Percentual ('%') or 'ppm' compostion (default is '%')
#     species : list of strings
#         List of chemical species.
#         Default is C(gr), N2, O2, H2, CO, CH4, CO2, H2O
    
#     Returns
#     -------
#     file : csv
#         Function return a CSV file as following data: %wt biomass ratio 
#         (assuming 1st fuel as coal), %wt moisture, T (degC), P (atm), 
#         equivalence ratio, steam-to-carbon ratio, O-to-C ratio, H-to-C ratio, 
#         species mole fractions, H2-to-CO ratio, % carbon conversion, 
#         gas yield (Nm3/kg), HHV (MJ/kg), % cold gas efficiency        
#     """
#     # convert all values to array
#     blend *= one
#     moisture *= one
#     T *= one
#     P *= one
#     air *= one
#     O2 *= one
#     ER *= one
#     steam *= one
#     SR *= one
#     # default values
#     steam_ = 0
#     SR_ = 0
#     air_ = 0
#     o2_ = 0
#     ER_ = 0
#     # get number of points
#     n_0 = np.size(fuel_id)
#     n_1 = np.size(blend)
#     n_2 = np.size(moisture)
#     n_3 = np.size(T)
#     n_4 = np.size(P)
    
#     if np.size(air) > 1:
#         n_5 = np.size(air)
#     elif np.size(O2) > 1:
#         n_5 = np.size(O2)
#     elif np.size(ER) > 1:
#         n_5 = np.size(ER)
#     else:
#         n_5 = 1
        
#     if np.size(steam) > 1:
#         n_6 = np.size(steam)
#     elif np.size(SR) > 1:
#         n_6 = np.size(SR)
#     else:
#         n_6 = 1
        
#     if format_ == 'ppm':
#         ft = 1e6
#     else:
#         ft = 1e2
# #    # start count minimum number of species
# #    minimum_species = []
#     # start calculations
#     for i in range(n_0-1): # asssumed 1st fuel as coal
#         csvfile = open(str(fuel_id[0]) + '-' + str(fuel_id[i+1]) + '.csv','w',
#                        newline='')
#         f = csv.writer(csvfile)
#         f.writerow(['% BR','% MC','T (C)','P (atm)','ER','SR','O/C',
#                     'H/C'] + species + ['H2/CO','% CC','Y (Nm3/kg)',
#                     'HHV (MJ/kg)','% CGE'])
#         for j in range(n_1): # %coal-biomass blend
#             frac = blending(blend[j], self[0,:], self[i+1,:])
#             for k in range(n_2): # moisture
#                 # get lhv to each moisture content of fuels
#                 fuel_lhv = feedstock.heating_values(fuel_id,moisture[k])['LHV']
#                 for l in range(n_3): # temperature
#                     for m in range(n_4): # pressure
#                         for o in range(n_5): # equivalence ratio
#                             if air.any() != 0:
#                                 air_ = air[o]
#                                 o2_ = 0
#                                 ER_ = equivalence_ratio(frac, 1.0, air[o])
#                             elif O2.any() != 0:
#                                 air_ = 0
#                                 o2_ = O2[o]
#                                 ER_ = equivalence_ratio(frac, 1.0, 0, O2[o])
#                             elif ER.any() != 0:
#                                 air_ = mass_of_air(frac, 1.0, ER[o])
#                                 o2_ = 0
#                                 ER_ = ER[o]
#                             for q in range(n_6): # steam-to-carbon ratio
#                                 if SR.any() != 0:
#                                     steam_ = mass_of_steam(frac, 1.0, SR[q])
#                                     SR_ = SR[q]
#                                 elif steam.any() != 0:
#                                     steam_ = steam[q]
#                                     SR_ = steam_to_carbon_ratio(frac, 1.0, 
#                                                                 steam[q])
#                                 hc,oc = ohc_ratio(frac, moisture[k], 1.0, 
#                                                   air_, o2_, steam_)
#                                 p,r = equilibrate_tp(frac, moisture[k], 1.0, 
#                                                      air_, o2_, steam_, 
#                                                      T[l]+273.15,
#                                                      ct.one_atm*P[m])
#                                 fuel_lhv_ = blending(blend[j], fuel_lhv[0], 
#                                                      fuel_lhv[i+1])
#                                 syngas_lhv_ = syngas_lhv(p, 1 + moisture[k])
#                                 eff = syngas_lhv_/fuel_lhv_
#                                 hhv = syngas_hhv(p, basis='fuel mass', 
#                                                  fuel_mass=1+moisture[k])
#                                 h2co = h2co_ratio(p)
#                                 cc = carbon_conversion(p,r)
#                                 y = gas_yield(p, basis='vol', db='y') # per kg
#                                 syngas = get_fraction(p, species, eps=small,
#                                                       db=db, 
#                                                       normalized=normalized)
#                                 f.writerow([100*blend[j], 100*moisture[k],
#                                             T[l], P[m], ER_,
#                                             SR_, oc, hc] + list(ft*syngas) 
#                                             + [h2co, 100*cc, y, hhv, 100*eff])
# #                                minimum_species = get_species(p, 
# #                                                              minimum_species, 
# #                                                              eps=1e-6)
#         csvfile.close()
# #        print minimum_species
#         print('Blend #' + str(i+1) + ' (' + str(fuel_id[0]) + '-' \
#               + str(fuel_id[i+1]) + '): DONE')
#     return None

# def coprocessing1(self, fuel_id, blend, moisture, T, P=1.0,
#                   air=0, O2=0, ER=0.4,
#                   steam=0, SR=0,
#                   small=None, db='n', normalized='n',
#                   species=['C(gr)','N2','O2','H2','CO','CH4','CO2','H2O']):
#     """
#     Cogasification calculations for binary blends of fuels.
    
#     Parameters
#     ----------
#     self : ndarray
#         Mass fraction of fuel #1 [kg/kg]
#     fuel_id : list of strings
#         List of ID fuel
#     blend : float|ndarray
#         Fuel #1 to fuel #2 ratio [kg/kg]
#     moisture : float|ndarray
#         Moisture mass fraction [kg/kg]
#     T : float|ndarray
#         Temperature [degC]
#     P : float|ndarray
#         Pressure [atm] (default is 1.0)
#     air : float|ndarray
#         Air amount [kg] (default is zero)
#     O2 : float|ndarray
#         O2 amount [kg] (default is zero)
#     ER : float|ndarray
#         Equivalence ratio [kmol/kmol]
#     steam : float|ndarray
#         Steam amount [kg] (default is zero)
#     SR : float|ndarray
#         Steam to carbon ratio [kmol/kmol] (default is zero)
#         basis: 1 kg coal-biomass blend, d.b.
#     small : float
#         Smallest number to report as a fraction value (default is None)
#     db : string
#         Get dry basis composition ('y') or web basis composition ('n') (default
#         is 'y')
#     normalized : string
#         Get normalized compostion ('y') or overall composition ('n') (default
#         is 'y')
#     species : list of strings
#         List of chemical species.
#         Default is N2, O2, H2, CO, CH4, CO2, H2O
    
#     Returns
#     -------
#     file : csv
#         Function return a CSV file as following data: %wt biomass ratio 
#         (assuming 1st fuel as coal), %wt moisture, T (degC), P (atm), 
#         equivalence ratio, steam-to-carbon ratio, O-to-C ratio, H-to-C ratio, 
#         species mole fractions, H2-to-CO ratio, % carbon conversion, 
#         gas yield (Nm3/kg), HHV (MJ/kg), % cold gas efficiency        
#     """
#     # convert all values to array
#     blend *= one
#     moisture = moisture*one
#     T *= one
#     air *= one
#     O2 *= one
#     ER *= one
#     steam *= one
#     SR *= one
#     # get number of points
#     n_0 = np.size(fuel_id)
#     n_1 = np.size(blend)    
#     # start calculations
#     for i in range(n_0-1): # asssumed 1st fuel as coal
#         csvfile = open(str(fuel_id[0]) + '-' + str(fuel_id[i+1]) + '.csv','w',
#                        newline='')
#         f = csv.writer(csvfile)
#         f.writerow(['% BR','% MC','T (C)','P (atm)','ER','SR','O/C',
#                     'H/C'] + species + ['H2/CO','% CC','Y (Nm3/kg)',
#                     'HHV (MJ/kg)','% CGE'])
#         for j in range(n_1): # %coal-biomass blend
#             frac = blending(blend[j], self[0,:], self[i+1,:])
#             # get lhv to each moisture content of fuels
#             fuel_lhv = feedstock.heating_values(fuel_id, moisture[j])['LHV']
#             if air.any() != 0:
#                 air_ = air[j]
#                 o2_ = 0
#                 ER_ = equivalence_ratio(frac, 1.0, air[j])
#             elif O2.any() != 0:
#                 air_ = 0
#                 o2_ = O2[j]
#                 ER_ = equivalence_ratio(frac, 1.0, 0, O2[j])
#             elif ER.any() != 0:
#                 air_ = mass_of_air(frac, 1.0, ER[j])
#                 o2_ = 0
#                 ER_ = ER[j]
#             else:
#                 air_ = 0
#                 o2_ = 0
#                 ER_ = 0
#             if SR.all() != 0:
#                 steam_ = mass_of_steam(frac, 1.0, SR[j])
#                 SR_ = SR[j]
#             elif steam.all() != 0:
#                 steam_ = steam[j]
#                 SR_ = steam_to_carbon_ratio(frac, 1.0, steam[j])
#             else:
#                 steam_ = 0
#                 SR_ = 0
#             hc, oc = ohc_ratio(frac, moisture[j], 1.0, air_, o2_, steam_)
#             p, r = equilibrate_tp(frac, moisture[j], 1.0, air_, o2_, steam_, 
#                                   T[j]+273.15, ct.one_atm*P)
#             fuel_lhv_ = blending(blend[j], fuel_lhv[0], fuel_lhv[i+1])
#             syngas_lhv_ = syngas_lhv(p, 1 + moisture[j])
#             eff = 100*syngas_lhv_/fuel_lhv_
#             hhv = syngas_hhv(p, basis='fuel mass', fuel_mass=1+moisture[j])
#             h2co = h2co_ratio(p)
#             cc = 100*carbon_conversion(p, r)
#             y = gas_yield(p, basis='vol', db='y') # per kg
#             syngas = get_fraction(p, species, eps=small, 
#                                   db=db, normalized=normalized)
#             f.writerow([100*blend[j], 100*moisture[j], T[j], P, ER_, SR_, 
#                         oc, hc] + list(100*syngas) + [h2co, cc, y, hhv, eff])
#         csvfile.close()
#         print('Blend #' + str(i+1) + ' (' + str(fuel_id[0]) + '-' \
#               + str(fuel_id[i+1]) + '): DONE')
#     return None

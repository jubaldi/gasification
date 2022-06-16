#!/usr/bin/env python

"""This script defines functions to equilibrium simulation of gasification 
processes. It uses some predefined functions from Cantera package.

@author = Rodolfo Rodrigues
@contact = rodolfo.rodrigues@ufsm.br
@data = April, 2012, rev.: June, 2013 (adapted to use cython Cantera)
"""

#==============================================================================
# import libraries/files
#==============================================================================
import pp as ppold
import pp2 as pp
import feedstock as fsold
import feedstock2 as fs
import fuel as fu
import outputs as op
import energy as en
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

def getFeed(fuelMix, moist=0.0, air=0.5, stm=0.0, airType='ER', stmType='SR'):
    '''
    This function creates a mixture of phases to denote the fuel.
    The fuel is composed as a mixture of char, gas, ash, and moisture phases.

    Parameters
    ----------
    fuelMix : Cantera 'Mixture' object
        Object containing the mole amount of each species in the dry fuel.
    moist : float
        Moisture content of fuel [kg/kg] in dry basis (default value is zero)
    air : float
        Mass amount of air [kg], equivalence ratio ER [kmol/kmol] or mass amount of pure O2 [kg] (default: ER=0.5)
    stm : float
        Mass amount of steam [kg] or steam to carbon ratio SR [kmol/kmol] (default: SR=0.0)
    airType : str
        Either 'ER', 'air' or 'O2' (default: 'ER')
    stmType : str
        Either 'SR' or 'steam' (default: 'SR')

    Returns
    -------
    feed : Cantera 'Mixture' object
        Object representing the mixture of phases in the feedstock.
    '''

    fuelMoles = fuelMix.species_moles # kmoles for each species

    totalFuelMass = fs.getFuelMass(fuelMix) # total fuel mass in kg

    moistMass = moist * totalFuelMass # mass of moisture in kg
    moistMoles = moistMass / pp.Mw['H2O'] # moles of moisture
    
    if stmType == 'SR':
        stmMass = fs.SRtosteam(fuelMix, stm)
        SR = stm
    elif stmType == 'steam':
        stmMass = stm
        SR = fs.steamtoSR(fuelMix, stm)
    else:
        raise ValueError('Invalid steam type')

    steamMoles = stmMass / pp.Mw['H2O'] # moles of steam
    H2OMoles = moistMoles + steamMoles # total moles of water
    
    if airType == 'ER':
        airMass = fs.ERtoair(fuelMix, air)
        pureO2Mass = 0
        ER = air
    elif airType == 'air':
        airMass = air
        pureO2Mass = 0
        ER = fs.airtoER(fuelMix, air)
    elif airType == 'O2':
        airMass = 0
        pureO2Mass = air
        ER = fs.airtoER(fuelMix, air=0.0, O2=air)
    else:
        raise ValueError('Invalid air type')

    airO2Moles = 0.23211606*airMass/pp.Mw['O2']
    airN2Moles = 0.75507754*airMass/pp.Mw['N2']
    airArMoles = 0.01280640*airMass/pp.Mw['Ar']

    pureO2Moles = pureO2Mass / pp.Mw['O2'] # moles of pure O2

    O2Moles = airO2Moles + pureO2Moles # total moles of O2

    # Creating feed mix
    feed = pp.mix()
    feedMoles = np.zeros(len(feed.species_names))
    # Adding species moles
    feedMoles += fuelMoles
    feedMoles[pp.i['H2O']] += H2OMoles
    feedMoles[pp.i['O2']] += O2Moles
    feedMoles[pp.i['N2']] += airN2Moles
    feedMoles[pp.i['Ar']] += airArMoles
    # Updating moles in mix
    feed.species_moles = feedMoles

    return feed

def isotGasification(fuelID, fuelMass=1.0, moist=0.0, T=1273.15, P=ct.one_atm, 
                    air=0.5, stm=0.0, airType='ER', stmType='SR', C_avail=1.0):
    '''
    Isothermal gasification calculation for a single fuel in a given condition.

    Parameters
    ----------
    fuelID : str
        ID of fuel as given by the database (fuels.csv)
    fuelMass : float
        The fuel mass [kg] (default: 1.0 kg)
    moist : float
        Moisture content of fuel [kg/kg] in dry basis (default value is zero)
    T : float
        Temperature [K] (default: 1273.15 K)
    P : float
        Pressure [Pa] (default: 101325 Pa)
    air : float
        Mass amount of air [kg], equivalence ratio ER [kmol/kmol] or mass amount of pure O2 [kg] (default: ER=0.5)
    stm : float
        Mass amount of steam [kg] or steam to carbon ratio SR [kmol/kmol] (default: SR=0.0)
    airType : str
        Either 'ER', 'air' or 'O2' (default: 'ER')
    stmType : str
        Either 'SR' or 'steam' (default: 'SR')
    C_avail : float
        Fraction of carbon from total that is available for equilibrium calculations (default: 1.0)

    Returns
    -------
    outlet : Cantera 'Mixture' object
        Object representing the mixture at equilibrium.
    inlet : Cantera 'Mixture' object
        Object representing the feed mixture.
    fuelMix : Cantera 'Mixture' object
        Object representing the dry fuel mixture.
    '''
    # Create fuel mix
    fuelMix = fs.getFuelMix(fuelID, fuelMass)

    # Separate out unavailable C moles
    fuelMoles = fuelMix.species_moles
    unavailable_C_moles = fuelMoles[pp.i["C(gr)"]] * (1 - C_avail)
    fuelMoles[pp.i["C(gr)"]] *= C_avail
    fuelMix.species_moles = fuelMoles

    # Create feed
    inlet = getFeed(fuelMix, moist, air, stm, airType, stmType)
    outlet = getFeed(fuelMix, moist, air, stm, airType, stmType)

    # Calculate equilibrium
    outlet.T = T
    outlet.P = P
    outlet.equilibrate('TP')

    # Add back unavailable C moles
    outMoles = outlet.species_moles
    outMoles[pp.i["C(gr)"]] += unavailable_C_moles
    outlet.species_moles = outMoles

    return outlet, inlet, fuelMix

def isotCogasification(fuel1, fuel2, mass=1.0, blend=0.5, moist=(0.0,0.0), T=1273.15, 
                    P=ct.one_atm, air=0.5, stm=0.0, airType='ER', stmType='SR'):
    '''
    Isothermal gasification calculation for a blend of 2 fuels in a given condition.

    Parameters
    ----------
    fuel1 : str
        ID of fuel #1 as given by the database (fuels.csv)
    fuel2 : str
        ID of fuel #2 as given by the database (fuels.csv)
    mass : float
        Total fuel blend mass [kg] (default: 1.0 kg)
    blend : float
        Fuel #2 to total fuel mass ratio [kg/kg] (default: 0.5)
    moist : tuple
        The moisture mass fraction for each fuel, in dry basis [kg/kg] (default: 0.0 for each fuel)
    T : float
        Temperature [K] (default: 1273.15 K)
    P : float
        Pressure [Pa] (default: 101325 Pa)
    air : float
        Mass amount of air [kg], equivalence ratio ER [kmol/kmol] or mass amount of pure O2 [kg] (default: ER=0.5)
    stm : float
        Mass amount of steam [kg] or steam to carbon ratio SR [kmol/kmol] (default: SR=0.0)
    airType : str
        Either 'ER', 'air' or 'O2' (default: 'ER')
    stmType : str
        Either 'SR' or 'steam' (default: 'SR')

    Returns
    -------
    outlet : Cantera 'Mixture' object
        Object representing the mixture at equilibrium.
    inlet : Cantera 'Mixture' object
        Object representing the feed mixture.
    fuelMix : Cantera 'Mixture' object
        Object representing the dry fuel mixture.
    '''
    # Create each fuel mixture
    fuel1Mass = mass*(1-blend)
    fuel2Mass = mass*blend
    fuel1Mix = fs.getFuelMix(fuel1, fuel1Mass)
    fuel2Mix = fs.getFuelMix(fuel2, fuel2Mass)

    # Create fuel blend
    fuelBlend = fs.blend(fuel1Mix, fuel2Mix)

    # Total moisture content [kg/kg]
    totalMoist = moist[0]*(1-blend) + moist[1]*blend

    # Create feed
    inlet = getFeed(fuelBlend, totalMoist, air, stm, airType='ER', stmType='SR')
    outlet = getFeed(fuelBlend, totalMoist, air, stm, airType='ER', stmType='SR')

    # Calculate equilibrium
    outlet.T = T
    outlet.P = P
    outlet.equilibrate('TP')

    return outlet, inlet, fuelBlend

def NonIsotGasification(fuelID, mass=1.0, moist=0.0, T0=1273.15, P=ct.one_atm,
                        air=0.5, stm=0.0, airType='ER', stmType='SR', heatLoss=0.0, guess=None):
    '''
    Non-isothermal gasification calculation for a single fuel in a given condition.
    Initial fuel temperature is always 25°C (298.15 K, T_ref).

    Parameters
    ----------
    fuelID : str
        ID of fuel as given by the database (fuels.csv)
    mass : float
        The fuel mass [kg] (default: 1.0 kg)
    moist : float
        Moisture content of fuel [kg/kg] in dry basis (default value is zero)
    T0 : float
        Initial air and steam temperature [K] (default: 1273.15 K)
    P : float
        Pressure [Pa] (default: 101325 Pa)
    air : float
        Mass amount of air [kg], equivalence ratio ER [kmol/kmol] or mass amount of pure O2 [kg] (default: ER=0.5)
    stm : float
        Mass amount of steam [kg] or steam to carbon ratio SR [kmol/kmol] (default: SR=0.0)
    airType : str
        Either 'ER', 'air' or 'O2' (default: 'ER')
    stmType : str
        Either 'SR' or 'steam' (default: 'SR')
    heatLoss : float
        Amount of heat lost to the environment [J]

    Returns
    -------
    outlet : Cantera 'Mixture' object
        Object representing the mixture at equilibrium.
    inlet : Cantera 'Mixture' object
        Object representing the feed mixture.
    fuelMix : Cantera 'Mixture' object
        Object representing the dry fuel mixture.
    '''
    if guess == None: guess = pp.To

    # Create fuel mix
    fuelMix = fs.getFuelMix(fuelID, mass)
    fuelMix.T = guess
    fuelMix.P = P

    # Create feed
    inlet = getFeed(fuelMix, moist, air, stm, airType, stmType)
    outlet = getFeed(fuelMix, moist, air, stm, airType, stmType)

    fuelMoles = sum(fuelMix.species_moles)
    HHV = fu.HV(fuelID, type='HHV', moist=moist)
    Hfo = en.hFormation(fuelID, HHV)

    moistMoles = moist*mass/pp.Mw['H2O']

    if stmType == 'SR':
        stmMass = fs.SRtosteam(fuelMix, stm)
        SR = stm
    elif stmType == 'steam':
        stmMass = stm
        SR = fs.steamtoSR(fuelMix, stm)
    else:
        raise ValueError('Invalid steam type')

    steamMoles = stmMass / pp.Mw['H2O'] # moles of steam
    
    if airType == 'ER':
        airMass = fs.ERtoair(fuelMix, air)
        pureO2Mass = 0
        ER = air
    elif airType == 'air':
        airMass = air
        pureO2Mass = 0
        ER = fs.airtoER(fuelMix, air)
    elif airType == 'O2':
        airMass = 0
        pureO2Mass = air
        ER = fs.airtoER(fuelMix, air=0.0, O2=air)
    else:
        raise ValueError('Invalid air type')

    airO2Moles = 0.23211606*airMass/pp.Mw['O2']
    airN2Moles = 0.75507754*airMass/pp.Mw['N2']
    airArMoles = 0.01280640*airMass/pp.Mw['Ar']

    pureO2Moles = pureO2Mass / pp.Mw['O2'] # moles of pure O2

    O2Moles = airO2Moles + pureO2Moles # total moles of O2
    totalMoles = fuelMoles + moistMoles + O2Moles + airN2Moles + airArMoles + steamMoles

    inlet_H = (fuelMoles*Hfo + moistMoles*(pp.Hfo['H2O(l)'] + pp.H_vap) +
                + O2Moles*pp.Hfo['O2'] + airN2Moles*pp.Hfo['N2'] + airArMoles*pp.Hfo['Ar'] +
                + steamMoles*pp.H_vap)/totalMoles # J/kmol
    
    desired_H = inlet_H*(1 - heatLoss) # J/kmol

    outlet.T = guess
    outlet.P = P
    outlet.equilibrate('TP')

    def residual(Temp):
        outlet.T = Temp
        outlet.P = P
        outlet_H = en.get_h_cp(outlet)
        return (outlet_H - desired_H)**2
    # estimate equilibrium temperature
    res = opt.minimize_scalar(residual,method='bounded',bounds=(200,6000),
                                bracket=(residual(1200),residual(3000)))
    # estimate equilibrium product composition
    Teq = res.x
    outlet.T = Teq
    outlet.P = P

    return outlet, inlet, fuelMix

def gasifier(fuelID, mass=1.0, moist=0.0, T=1273.15, P=ct.one_atm, 
                air=0.5, stm=0.0, airType='ER', stmType='SR', C_avail=1.0, isot=True,
                species=['C(gr)','N2','O2','H2','CO','CH4','CO2','H2O']):
    '''
    Creates a full report of the outputs for a given gasification condition.

    Parameters
    ----------
    fuelID : str
        ID of fuel as given by the database (fuels.csv)
    mass : float
        The fuel mass [kg] (default: 1.0 kg)
    moist : float
        Moisture content of fuel [kg/kg] in dry basis (default value is zero)
    T : float
        Temperature [K] (default: 1273.15 K)
    P : float
        Pressure [Pa] (default: 101325 Pa)
    air : float
        Mass amount of air [kg], equivalence ratio ER [kmol/kmol] or mass amount of pure O2 [kg] (default: ER=0.5)
    stm : float
        Mass amount of steam [kg] or steam to carbon ratio SR [kmol/kmol] (default: SR=0.0)
    airType : str
        Either 'ER', 'air' or 'O2' (default: 'ER')
    stmType : str
        Either 'SR' or 'steam' (default: 'SR')
    C_avail : float
        Fraction of carbon from total that is available for equilibrium calculations (default: 1.0)
    isot : bool
        If True, use isothermal gasification calculation.
    species : list
        List of species to be included in the report.
        (default: ['C(gr)','N2','O2','H2','CO','CH4','CO2','H2O'])

    Returns
    -------
    report : dict
        Dictionary containing the report.
    '''
    if isot:
        # Isothermal gasification
        outlet, inlet, fuelMix = isotGasification(fuelID, mass, moist, T, P, air, stm, airType, stmType, C_avail)
    else:
        raise NotImplementedError('Non-isothermal gasification not yet implemented.')

    if stmType == 'SR':
        stmMass = fs.SRtosteam(fuelMix, stm)
        SR = stm
    elif stmType == 'steam':
        stmMass = stm
        SR = fs.steamtoSR(fuelMix, stm)
    else:
        raise ValueError('Invalid steam type')
    
    if airType == 'ER':
        airMass = fs.ERtoair(fuelMix, air)
        pureO2Mass = 0
        ER = air
    elif airType == 'air':
        airMass = air
        pureO2Mass = 0
        ER = fs.airtoER(fuelMix, air)
    elif airType == 'O2':
        airMass = 0
        pureO2Mass = air
        ER = fs.airtoER(fuelMix, air=0.0, O2=air)
    else:
        raise ValueError('Invalid air type')

    # Create report
    report = {}

    report['FuelID'] = fuelID
    report['Fuel'] = fu.fuels.loc[fuelID]['Description']
    report['Fuel mass'] = mass
    report['Moisture'] = moist
    report['T'] = T - 273.15
    report['P'] = P/ct.one_atm
    report['ER'] = ER
    report['SR'] = SR

    OC, HC = fs.OHCratio(fuelMix)
    report['O/C'] = OC
    report['H/C'] = HC

    fracs = op.getAmounts(outlet, species, norm=True)
    for i, s in enumerate(species):
        report[s] = fracs[i]

    report['H2/CO'] = op.H2CO(outlet)
    report['CC'] = op.carbonConv(outlet, inlet)*100
    report['Y'] = op.gasYield(outlet, basis='vol')/mass
    report['HHV'] = op.syngasHHV(outlet, basis='fuel mass', fuelMass=mass)

    fuelLHV = fu.HV(fuelID, type='LHV', moist=moist)
    report['CGE'] = op.coldGasEff(outlet, fuelLHV, moist=moist)*100  

    return report

def cogasifier(fuel1, fuel2, mass=1.0, blend=0.5, moist=(0.0,0.0), T=1273.15, 
                P=ct.one_atm, air=0.5, stm=0.0, airType='ER', stmType='SR', isot=True,
                species=['C(gr)','N2','O2','H2','CO','CH4','CO2','H2O']):
    '''
    Creates a full report of the outputs for a given cogasification condition.

    Parameters
    ----------
    fuel1 : str
        ID of fuel #1 as given by the database (fuels.csv)
    fuel2 : str
        ID of fuel #2 as given by the database (fuels.csv)
    mass : float
        Total fuel blend mass [kg] (default: 1.0 kg)
    blend : float
        Fuel #2 to total fuel mass ratio [kg/kg] (default: 0.5)
    moist : tuple
        The moisture mass fraction for each fuel, in dry basis [kg/kg] (default: 0.0 for each fuel)
    T : float
        Temperature [K] (default: 1273.15 K)
    P : float
        Pressure [Pa] (default: 101325 Pa)
    air : float
        Mass amount of air [kg], equivalence ratio ER [kmol/kmol] or mass amount of pure O2 [kg] (default: ER=0.5)
    stm : float
        Mass amount of steam [kg] or steam to carbon ratio SR [kmol/kmol] (default: SR=0.0)
    airType : str
        Either 'ER', 'air' or 'O2' (default: 'ER')
    stmType : str
        Either 'SR' or 'steam' (default: 'SR')
    isot : bool
        If True, use isothermal gasification calculation.
    species : list
        List of species to be included in the report.
        (default: ['C(gr)','N2','O2','H2','CO','CH4','CO2','H2O'])

    Returns
    -------
    report : dict
        Dictionary containing the report.
    '''
    if isot:
        # Isothermal gasification
        outlet, inlet, fuelBlend = isotCogasification(fuel1, fuel2, mass, blend, moist, T, P, air, stm, airType, stmType)
    else:
        raise NotImplementedError('Non-isothermal gasification not yet implemented.')

    if stmType == 'SR':
        stmMass = fs.SRtosteam(fuelBlend, stm)
        SR = stm
    elif stmType == 'steam':
        stmMass = stm
        SR = fs.steamtoSR(fuelBlend, stm)
    else:
        raise ValueError('Invalid steam type')
    
    if airType == 'ER':
        airMass = fs.ERtoair(fuelBlend, air)
        pureO2Mass = 0
        ER = air
    elif airType == 'air':
        airMass = air
        pureO2Mass = 0
        ER = fs.airtoER(fuelBlend, air)
    elif airType == 'O2':
        airMass = 0
        pureO2Mass = air
        ER = fs.airtoER(fuelBlend, air=0.0, O2=air)
    else:
        raise ValueError('Invalid air type')

    # Create report
    report = {}

    report['Fuel #1'] = fu.fuels.loc[fuel1]['Description']
    report['Fuel #2'] = fu.fuels.loc[fuel2]['Description']
    report['Mass'] = mass
    report['Blend'] = blend*100

    blendMoist = moist[0]*(1-blend) + moist[1]*blend
    report['Moisture'] = blendMoist

    report['T'] = T - 273.15
    report['P'] = P/ct.one_atm
    report['ER'] = ER
    report['SR'] = SR

    OC, HC = fs.OHCratio(inlet)
    report['O/C'] = OC
    report['H/C'] = HC

    fracs = op.getAmounts(outlet, species, norm=True)
    for i, s in enumerate(species):
        report[s] = fracs[i]*100

    report['H2/CO'] = op.H2CO(outlet)
    report['CC'] = op.carbonConv(outlet, inlet)*100
    report['Y'] = op.gasYield(outlet, basis='vol')/mass
    report['HHV'] = op.syngasHHV(outlet, basis='fuel mass', fuelMass=mass)

    fuel1LHV = fu.HV(fuel1, type='LHV', moist=moist[0])
    fuel2LHV = fu.HV(fuel2, type='LHV', moist=moist[1])
    blendLHV = fuel1LHV*(1-blend) + fuel2LHV*blend
    report['CGE'] = op.coldGasEff(outlet, blendLHV, moist=blendMoist)*100 

    return report

def findTquasi(fuelID, experimental, mass=1.0, moist=0.0, T0=1273.15, P=ct.one_atm, 
                air=0.5, stm=0.0, airType='ER', stmType='SR',
                species=['C(gr)','N2','O2','H2','CO','CH4','CO2','H2O']):
    '''
    Finds the quasi-equilibrium temperature for a given gasification condition.
    The quasi-equilibrium temperature is the temperature that minimizes the square error of compositions
    when applied to the equilibrium model.
    Must be supplied with an array of experimental compositions.

    Parameters
    ----------
    fuelID : str
        ID of fuel as given by the database (fuels.csv)
    experimental : array
        Array of measured experimental gas compositions (mole fractions). Must be of same length as species.
    mass : float
        The fuel mass [kg] (default: 1.0 kg)
    moist : float
        Moisture content of fuel [kg/kg] in dry basis (default value is zero)
    T0 : float
        Temperature [K], will be used as first guess (default: 1273.15 K)
    P : float
        Pressure [Pa] (default: 101325 Pa)
    air : float
        Mass amount of air [kg], equivalence ratio ER [kmol/kmol] or mass amount of pure O2 [kg] (default: ER=0.5)
    stm : float
        Mass amount of steam [kg] or steam to carbon ratio SR [kmol/kmol] (default: SR=0.0)
    airType : str
        Either 'ER', 'air' or 'O2' (default: 'ER')
    stmType : str
        Either 'SR' or 'steam' (default: 'SR')
    species : list
        List of species to be used for error minimization.
        (default: ['C(gr)','N2','O2','H2','CO','CH4','CO2','H2O'])

    Returns
    -------
    Tquasi : float
        Quasi-equilibrium temperature [K]
    sqerr : float
        Sum of square errors.
    '''
    # Defining objective function
    def error(T):
        err = 0
        predicted = np.zeros(len(species))
        report = gasifier(fuelID, mass=mass, moist=moist, T=T, P=P, air=air, 
                 stm=stm, airType=airType, stmType=stmType, isot=True, species=species)
        for i, sp in enumerate(species):
            predicted[i] = report[sp]
            err += (predicted[i] - experimental[i])**2
        return err
    
    # Applying minimize
    Tquasi = opt.minimize(error, T0).x[0]
    sqerr = error(Tquasi)
    return Tquasi, sqerr

def findC_unav(fuelID, experimental, mass=1.0, moist=0.0, T0=1273.15, P=ct.one_atm, 
                air=0.5, stm=0.0, airType='ER', stmType='SR',
                species=['C(gr)','N2','O2','H2','CO','CH4','CO2','H2O']):
    def error(C_unav):
        err = 0
        predicted = np.zeros(len(species))
        report = gasifier(fuelID, mass=mass, moist=moist, T=T, P=P, air=air, 
                 stm=stm, airType=airType, stmType=stmType, isot=True, species=species)
        for i, sp in enumerate(species):
            predicted[i] = report[sp]
            err += (predicted[i] - experimental[i])**2
        return err
    
    # Applying minimize
    Tquasi = opt.minimize(error, T0).x[0]
    sqerr = error(Tquasi)
    return Tquasi, sqerr

# def oneToOne(conditions):
#     '''
#     Given a list of multiple items, some of them lists, returns a list of conditions matching 1 to 1 each item.

#     Parameters
#     ----------
#     conditions : list
#         List of conditions, where each item can be a list or a single value. All lists must have the same length.
    
#     Returns
#     -------
#     oneToOne : list
#         List of conditions matching 1 to 1 each item.
#     '''

#     # Find first item that is of type list
#     for c in conditions:
#         if isinstance(c, list):
#             Len = len(c)
    
#     listIndex = []
#     # Now, check if all list items are of the same length
#     for c in conditions:
#         if isinstance(c, list):
#             if len(c) != Len:
#                 raise ValueError('All lists must have the same length.')

#     oneToOne = []
#     for j in range(Len):
#         item = []
#         for c in conditions:
#             if isinstance(c, list):
#                 item.append(c[j])
#             else:
#                 item.append(c)
#         oneToOne.append(item)

#     return oneToOne
        
# def combinatorial(conditions):
#     '''
#     Given a list of multiple items, some of them lists, returns a list of conditions combining all of them.

#     Parameters
#     ----------
#     conditions : list
#         List of conditions, where each item can be a list or a single value.
    
#     Returns
#     -------
#     combine : list
#         List of conditions combining all items.
#     '''
#     combine = []
#     for c in conditions:
#         if isinstance(c, list):
#             item = []
#             for i, cond in enumerate(c):
#                 combine.append(item)
#             item.append(c[j])
#         else:
#             item.append(c)
#     oneToOne.append(item)

#     return None
            
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

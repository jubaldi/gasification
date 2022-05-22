#!/usr/bin/env python

"""This script defines functions to create inputs for the gasification simulations.

@author = Rodolfo Rodrigues
@contact = rodolfo.rodrigues@ufsm.br
@data = April, 2012, rev.: June, 2013 (adapted to use cython Cantera)
"""

# FIXME: 2013/05/25 Not possible to use a list of single value as parameter to 
# functions.
# QSTN: 2022/04/15 Ask professor if no data should return 0 or np.nan.
# TODO: 2022/04/27 Create functions: from fuelID and mass get Mix, from Mix get mass.
# TODO: 2022/04/18 Biochemical composition should be in csv.
# TODO: 2022/05/01 Find ash composition for petroleum coke.
# TODO: 2022/05/01 Add a cantera mixture for steam + water.
# TODO: 2022/05/01 Find heating values HHV and LHV from correlations.
# QSTN: 2022/05/01 Should we assume that remaining mass is fixed carbon?

#==============================================================================
# import libraries
#==============================================================================
from matplotlib.pyplot import get
import numpy as np
import pandas as pd
import pp2 as pp
import fuel as fu

def getFuelMix(fuelID, fuelMass):
    '''
    Create a Cantera 'Mixture' object representing the given fuel.

    Parameters
    ----------
    fuelID : str
        The fuel ID as given by the CSV.
    fuelMass : float
        The fuel mass [kg]

    Returns
    -------
    fuelMix : Cantera 'Mixture' object
        The fuel mixture object.
    '''
    comp = fu.fuelComp(fuelID)
    massBySpecies = {key: value*fuelMass for key, value in comp.items()} # kg

    molesBySpecies = {key: value/pp.Mw[key] for key, value in massBySpecies.items()} # kmol

    fuelMix = pp.mix()

    molesInOrder = np.zeros(len(fuelMix.species_names))
    for index, species in enumerate(fuelMix.species_names):
        if species in molesBySpecies.keys():
            molesInOrder[index] = molesBySpecies[species]
        else:
            molesInOrder[index] = 0

    fuelMix.species_moles = molesInOrder

    return fuelMix


def getFuelMass(fuelMix):
    '''
    Calculate the mass of the given fuel mixture.

    Parameters
    ----------
    fuelMix : Cantera 'Mixture' object
        The fuel mixture object.

    Returns
    -------
    fuelMass : float
        The fuel mass [kg]
    '''
    fuelMass = 0
    for index, species in enumerate(fuelMix.species_names):
        fuelMass += fuelMix.species_moles[index]*pp.Mw[species] # kg

    return fuelMass

def stoichO2(fuelMix):
    '''
    Calculates the stoichometric amount of oxygen moles required to
    combust the fuel completely. This is required to calculate the
    equivalence ratio (ER).

    Parameters
    ----------
    fuelMix : Cantera 'Mixture' object
        The fuel mixture object.

    Returns
    -------
    stoichO2frac : float
        The stoichometric amount of O2 in moles by the fuel moles [kmol/kmol]
    '''

    totalMoles = sum(fuelMix.species_moles)

    m = lambda e: fuelMix.element_moles(e)
    stoicO2Mole = (m('C') + 0.25*m('H') - 0.5*m('O') + m('S') - 0.5*m('Cl') +
                    + m('Si') + 0.5*m('Ca') + 1.5*m('Al') +
                    + 1.5*m('Fe') + 0.25*m('Na') + 0.25*m('K') + 0.5*m('Mg') +
                    + 2.5*m('P') + m('Ti') + 1.5*m('Cr'))

    stoichO2frac = stoicO2Mole / totalMoles

    return stoichO2frac

def airtoER(fuelMix, air, O2=0.0):
    '''
    Calculates the equivalence ratio (ER) of the fuel mixture for the
    given air / pure O2 mass.

    Parameters
    ----------
    fuelMix : Cantera 'Mixture' object
        Fuel mixture object.
    air : float
        Air mass [kg]
    O2 : float
        Pure O2 mass [kg]

    Returns
    -------
    ER : float
        Equivalence ratio [kmol/kmol]
    '''
    totalMoles = sum(fuelMix.species_moles)
    stoich = stoichO2(fuelMix)
    stoichO2Moles = stoich * totalMoles

    airMoles = air / pp.Mw_air
    pureO2Moles = O2 / pp.Mw['O2']

    airO2Moles = 0.21 * airMoles

    ER = (airO2Moles + pureO2Moles) / stoichO2Moles

    return ER

def ERtoair(fuelMix, ER=1.0):
    '''
    Calculates the air mass required to combust the fuel mixture for the
    given equivalence ratio.
    
    Parameters
    ----------
    fuelMix : Cantera 'Mixture' object
        Fuel mixture object.
    ER : float
        Equivalence ratio [kmol/kmol]
    
    Returns
    -------
    air : float
        Air mass [kg]
    '''
    totalMoles = sum(fuelMix.species_moles)
    stoich = stoichO2(fuelMix)
    stoichO2Moles = stoich * totalMoles

    realO2Moles = ER * stoichO2Moles

    airMoles = realO2Moles / 0.21

    air = airMoles * pp.Mw_air

    return air

def steamtoSR(mix, steam):
    '''
    Calculates the Steam to Carbon Ratio (SR) equivalent to the given steam mass for the mixture.

    Parameters
    ----------
    mix : Cantera 'Mixture' object
        Fuel mixture object.
    steam : float
        Steam mass [kg]

    Returns
    -------
    SR : float
        Steam to Carbon Ratio [kmol/kmol]
    '''
    steamMoles = steam / pp.Mw['H2O']
    if mix.element_moles('C') == 0:
        raise ValueError('No carbon in mixture.')
    SR = steamMoles / mix.element_moles('C')
    return SR

def SRtosteam(mix, SR):
    '''
    Calculates the steam mass equivalent to the given Steam to Carbon Ratio (SR) for the mixture.

    Parameters
    ----------
    mix : Cantera 'Mixture' object
        Fuel mixture object.
    SR : float
        Steam to Carbon Ratio [mol/mol]
    
    Returns
    -------
    steam : float
        Steam mass [kg]
    '''
    if mix.element_moles('C') == 0:
        raise ValueError('No carbon in mixture.')
    steamMoles = SR * mix.element_moles('C')
    steam = steamMoles * pp.Mw['H2O']
    return steam

def OHCratio(mix):
    '''
    Calculates H/C and O/C ratios for given mixture.

    Parameters
    ----------
    mix : Cantera 'Mixture' object
        Fuel mixture object.

    Returns
    -------
    OC : float
        O/C ratio [kmol/kmol]
    HC : float
        H/C ratio [kmol/kmol]
    '''
    if mix.element_moles('C') == 0:
        raise ValueError('No carbon in mixture.')
    HC = mix.element_moles('H') / mix.element_moles('C')
    OC = mix.element_moles('O') / mix.element_moles('C')
    return OC, HC

def blend(fuelMix1, fuelMix2):
    '''
    Blend two fuel mixtures, using the mass already given by Cantera.

    Parameters
    ----------
    fuelMix1 : Cantera 'Mixture' object
        Fuel mixture object.
    fuelMix2 : Cantera 'Mixture' object
        Fuel mixture object.
    
    Returns
    -------
    fuelBlend : Cantera 'Mixture' object
        Mixture object representing fuel blend.
    '''
    moles1 = fuelMix1.species_moles
    moles2 = fuelMix2.species_moles
    moles = moles1 + moles2
    fuelBlend = pp.mix()
    fuelBlend.species_moles = moles
    return fuelBlend

# def mass_to_mole_fraction(self, Mw1, Mw2):
#     """
#     Convert mass fraction to mole fraction for dual-fuel blends.
    
#     Parameters
#     ----------
#     self : ndarray
#         Mass fraction of fuel #1 [kg/kg]
#     Mw1 : float
#         Molecular weight of fuel #1 [kg/kmol]
#     Mw2 : float
#         Molecular weight of fuel #2 [kg/kmol]
    
#     Returns
#     -------
#     mole_fraction : ndarray
#         Mole fraction of fuel #1 [kmol/kmol]
#     """
#     idx = len(self)
#     if (self.ndim == 1):
#         mole_fraction = self/Mw1/(self/Mw1 + (1.0 - self)/Mw2)
#     else:
#         mole_fraction = np.zeros(idx,'d')
#         for i in range(idx):
#             mole_fraction[i] = self[i]/Mw1/(self[i]/Mw1 + (1 - self[i])/Mw2)
#     return mole_fraction

# def mole_to_mass_fraction(self, Mw1, Mw2):
#     """
#     Convert mole fraction to mass fraction for dual-fuel blends.
    
#     Parameters
#     ----------
#     self : ndarray
#         Mole fraction of fuel #1 [kmol/kmol]
#     Mw1 : float
#         Molecular weight of fuel #1 [kg/kmol]
#     Mw2 : float
#         Molecular weight of fuel #2 [kg/kmol]
    
#     Returns
#     -------
#     mass_fraction : ndarray
#         Mass fraction of fuel #1 [kg/kg]
#     """
#     idx = len(self)
#     if (self.ndim == 1):
#         mass_fraction = self*Mw1/(Mw2 - self(Mw1 - Mw2))
#     else:
#         mass_fraction = np.zeros(idx,'d')
#         for i in range(idx):
#             mass_fraction[i] = self[i]*Mw1/(Mw2 - self[i]*(Mw1 - Mw2))
#     return mass_fraction

# def mixture(f, prop1, prop2):
#     n1 = np.size(f)
#     if (prop1.ndim <= 0):
#         prop3 = np.zeros((n1))
#         for i in range(n1):
#             prop3[i] = f[i]*prop1 + (1.0 - f[i])*prop2
#     else:
#         n2 = len(prop1)
#         prop3 = np.zeros((n1,n2))
#         for i in range(n1):
#             for j in range(n2):
#                 prop3[i,j] = f[i]*prop1[j] + (1.0 - f[i])*prop2[j]
#     return prop3

# def blending(f, coal, biomass):
#     """
#     f : float
#         %wt biomass in coal-biomass blend
#     """
#     return (1.0 - f)*coal + (f)*biomass

# def avg_error(mes, sim):
#     """
#     Return average error
#     sim : ndarray
#         simulated values
#     mes: ndarray
#         mesuared values
#     """
#     return np.sum(np.abs(sim-mes)/mes)/len(mes)

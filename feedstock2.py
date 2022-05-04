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

# TODO LATER: 2022/05/01 Get heat of formation.

#==============================================================================
# import libraries
#==============================================================================
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
    massBySpecies = {key: value*fuelMass for key, value in comp.items()}

    molesBySpecies = {key: value*1000/pp.Mw[key] for key, value in massBySpecies.items()}

    fuelMix = pp.f

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
        fuelMass += fuelMix.species_moles[index]*pp.Mw[species]/1000

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
        The stoichometric amount of O2 in moles by the fuel moles [mol/mol]
    '''

    totalMoles = sum(fuelMix.species_moles)

    m = lambda e: fuelMix.element_moles(e)
    stoicO2Mole = m('C') + 0.25*m('H') - 0.5*m('O') + m('S') \
                - 0.5*m('Cl') + m('Si') + 0.5*m('Ca') + 1.5*m('Al') \
                + 1.5*m('Fe') + 0.25*m('Na') + 0.25*m('K') + 0.5*m('Mg') \
                + 2.5*m('P') + m('Ti') + 1.5*m('Cr') \

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
        Equivalence ratio [dimensionless]
    '''
    totalMoles = sum(fuelMix.species_moles)
    stoich = stoichO2(fuelMix)
    stoichO2Moles = stoich * totalMoles

    airMoles = air*1000 / pp.Mw_air
    pureO2Moles = O2*1000 / pp.Mw['O2']

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
        Equivalence ratio [dimensionless]
    
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

    air = airMoles * pp.Mw_air / 1000

    return air

# def steam_to_carbon_ratio(self, fuel, steam):
#     mol = chon_moles(self, 0, fuel, 0, 0, 0)
#     mol_of_C = mol[0]
#     mol_of_steam = steam / pp.Mw[pp.i_H2O]
#     return mol_of_steam / mol_of_C
    
# def mass_of_steam(self, fuel, SR=0):
#     mol = chon_moles(self, 0, fuel, 0, 0, 0)
#     mol_of_C = mol[0]
#     mol_of_steam = SR * mol_of_C
#     return mol_of_steam * pp.Mw[pp.i_H2O]
    
# def ohc_ratio(self, moist, fuel, air, o2, stm):
#     C, H, O, N = chon_moles(self, moist, fuel, air, o2, stm)
#     return H/C, O/C

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

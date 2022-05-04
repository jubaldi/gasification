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
import fuels as fu

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

# def get_fuel_db(self):
# #    fuel = get_feed(self, zero, one, zero) # 1 kg of fuel in d.b.
#     fuel = get_feed(self) # 1 kg of fuel in d.b.
#     nsp = fuel.n_species
#     sp = fuel.species_moles
#     # initiate variables
#     mol_of_C = 0
#     mol_of_H = 0
#     mol_of_O = 0
#     mol_of_S = 0
#     mol_of_Cl = 0
#     mol_of_Si = 0
#     mol_of_Ca = 0
#     mol_of_Al = 0
#     mol_of_Fe = 0
#     mol_of_Na = 0
#     mol_of_K = 0
#     mol_of_Mg = 0
#     mol_of_P = 0
#     mol_of_Ti = 0
#     mol_of_Cr = 0
#     mol_of_Ar = 0
#     mol = 0
#     # count moles of C,H,O in fuel species
#     # IMPORTANT: I have to count S, Cl and ash species for precise estimation 
#     # of stoichiometric oxygen amount. This is important mainly for high ash
#     # fuels
#     for i in range(nsp):
#         if sp[i] != 0:
# #            if i != fuel.species_index('gas', 'H2O'):
# #                if i != fuel.species_index('gas', 'CO2'):
# #                    mol_of_C += sp[i] * fuel.n_atoms(i, 'C')
# #                    mol_of_H += sp[i] * fuel.n_atoms(i, 'H')
# #                    mol_of_O += sp[i] * fuel.n_atoms(i, 'O')
#             mol_of_C += sp[i] * fuel.n_atoms(i, 'C')
#             mol_of_H += sp[i] * fuel.n_atoms(i, 'H')
#             mol_of_O += sp[i] * fuel.n_atoms(i, 'O')
#             mol_of_S += sp[i] * fuel.n_atoms(i, 'S')
#             mol_of_Cl += sp[i] * fuel.n_atoms(i, 'Cl')
#             mol_of_Si += sp[i] * fuel.n_atoms(i, 'Si')
#             mol_of_Ca += sp[i] * fuel.n_atoms(i, 'Ca')
#             mol_of_Al += sp[i] * fuel.n_atoms(i, 'Al')
#             mol_of_Fe += sp[i] * fuel.n_atoms(i, 'Fe')
#             mol_of_Na += sp[i] * fuel.n_atoms(i, 'Na')
#             mol_of_K += sp[i] * fuel.n_atoms(i, 'K')
#             mol_of_Mg += sp[i] * fuel.n_atoms(i, 'Mg')
#             mol_of_P += sp[i] * fuel.n_atoms(i, 'P')
#             mol_of_Ti += sp[i] * fuel.n_atoms(i, 'Ti')
#             mol_of_Cr += sp[i] * fuel.n_atoms(i, 'Cr')
#             mol_of_Ar += sp[i] * fuel.n_atoms(i, 'Ar')
#             mol += sp[i]
#     # normalise per mole of fuel
#     mol_of_C /= mol
#     mol_of_H /= mol
#     mol_of_O /= mol
#     mol_of_S /= mol
#     mol_of_Cl /= mol
#     mol_of_Si /= mol
#     mol_of_Ca /= mol
#     mol_of_Al /= mol
#     mol_of_Fe /= mol
#     mol_of_Na /= mol
#     mol_of_K /= mol
#     mol_of_Mg /= mol
#     mol_of_P /= mol
#     mol_of_Ti /= mol
#     mol_of_Cr /= mol
#     mol_of_Ar /= mol
#     # stoichiometric moles of oxygen per mole of fuel
#     stoic = mol_of_C + 0.25*mol_of_H - 0.5*mol_of_O + mol_of_S \
#             - 0.5*mol_of_Cl + mol_of_Si + 0.5*mol_of_Ca + 3/2*mol_of_Al \
#             + 3/2*mol_of_Fe + 0.25*mol_of_Na + 0.25*mol_of_K + 0.5*mol_of_Mg \
#             + 2.5*mol_of_P + mol_of_Ti + 3/2*mol_of_Cr
#     if stoic < 0:   # FIXME: Figure out the issue of a negative stoic
#                     # oxygen. This happens when there is a fuel with high
#                     # oxygen content, that is, 
#                     # 0.5*mol_of_O > mol_of_C + 0.25*mol_of_H
#         stoic += 0.5*mol_of_O
#     return fuel, stoic

# def mass_of_air(self, fuel, ER=1.0):
#     fuel_db, stoic = get_fuel_db(self)
#     mol_of_fuel = fuel * np.sum(self/pp.Mw_f[:-1])
#     # mole amount of gasifying agent
#     mol_of_air = ER * stoic * mol_of_fuel/0.21
#     # mass amount of gasifying agent
#     return mol_of_air * pp.Mw_air

# def equivalence_ratio(self, fuel, air, o2=0):
#     fuel_db, stoic = get_fuel_db(self)
#     mol_of_fuel = fuel * np.sum(self/pp.Mw_f[:-1])
#     if air!=0 and o2==0:
#         mol_of_O2 = 0.21 * (air/pp.Mw_air)
#     elif air==0 and o2!=0:
#         mol_of_O2 = o2/pp.Mw[pp.i_O2]
#     else:
#         mol_of_O2 = 0.21 * (air/pp.Mw_air) + o2/pp.Mw[pp.i_O2]
#     return mol_of_O2/(stoic * mol_of_fuel)

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

# def chon_moles(self, moist, fuel, air, o2, stm):
#     f = get_feed(self, moist, fuel, air, o2, stm)
#     nsp = f.n_species
#     sp = f.species_moles
#     # initiate variables
#     mol_of_C = 0
#     mol_of_H = 0
#     mol_of_O = 0
#     mol_of_N = 0
#     # count moles of C,H,O in fuel species
#     for i in range(nsp):
#         if sp[i] != 0:
#             mol_of_C += sp[i] * f.n_atoms(i, 'C')
#             mol_of_H += sp[i] * f.n_atoms(i, 'H')
#             mol_of_O += sp[i] * f.n_atoms(i, 'O')
#             mol_of_N += sp[i] * f.n_atoms(i, 'N')
#     return mol_of_C, mol_of_H, mol_of_O, mol_of_N
    
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

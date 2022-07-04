#!/usr/bin/env python

"""
This script defines functions to process outputs after
equilibrium simulation of gasification processes are done. 
It uses some predefined functions from Cantera package.

@author = Rodolfo Rodrigues
@contact = rodolfo.rodrigues@ufsm.br
@data = May, 2022

"""

#==============================================================================
# import libraries/files
#==============================================================================
import sys
import cantera as ct
import numpy as np

import pp
import feedstock as fs
import fuel as fu
import gasifier as g
import energy as en
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

def gasYield(outlet, basis='vol', db=True):
    """
    Gas yield of reactor outlet.

    Parameters
    ----------
    outlet : Cantera 'Mixture' object
        Gasifier outlet object.
    basis : string
        'mole' = mole amount in kmol
        'mass' = mass amount in kg
        'vol' = normal volume amount, Nm³
        Normal condition at 273.15K and 1 atm.
    db : string
        Dry basis ('y', default) or wet basis ('n')
    
    Returns
    -------
    yield : float
        Syngas yield [kmol] [kg] [Nm³]
    """
    # mole of gas species
    moles = outlet.phase_moles(p='gas')

    mass = 0
    for i, sp in enumerate(outlet.phase(outlet.phase_index('gas')).species_names):
        mass += outlet.species_moles[pp.i[sp]]*pp.Mw[sp]

    if db:
        moles -= outlet.species_moles[pp.i['H2O']]
        mass -= outlet.species_moles[pp.i['H2O']] * pp.Mw['H2O']

    moles -= outlet.species_moles[pp.i['N2']]
    mass -= outlet.species_moles[pp.i['N2']]*pp.Mw['N2']
    vol = moles*R*Tn/Pn
     
    if basis == 'mole':
        return moles
    elif basis == 'mass':
        return mass
    elif basis == 'vol':
        return vol
    else:
        raise ValueError('Basis must be mole, mass or vol')

def getSpecies(mix, species=[], eps=1e-6):
    '''
    Gets a list of species which mole fractions in 'mix' are higher than 'eps'.
    this function is useful to find a minimum number of species to handle out a
    chemical equilibrium problem.
    '''
    i = 1
    fracs = mix.species_moles / sum(mix.species_moles)

    while i < pp.nsp:
        if fracs[i] > eps:
            speciesName = pp.f.species_name(i)
            if ('L)' not in speciesName) and (speciesName not in species):
                species.append(speciesName)
        i += 1
    return species

def getAmounts(mix, species, norm=False, db=False, eps=None, mass=False, phase='gas'):
    '''
    Gets an array of amounts (mole, mass, or fractions) of the given species.

    Parameters
    ----------
    mix : Cantera 'Mixture' object
        Fuel mixture object.
    species : list
        List of species names.
    norm : boolean
        If True, returns normalized amounts, fractions. (default: False, absolute amounts)
    db : boolean
        If True, returns amounts in dry basis. (default: False, wet basis)
    eps : float
        Any values smaller than 'eps' become zero. (default: None)
    mass : boolean
        If True, returns mass amounts [kg]. (default: False, mole amounts [kmol])
    phase : str
        Phase of the mixture, for fraction calculation as total amount: 'gas', 'solid', or all. (default: 'gas')

    Returns
    -------
    amounts : array
        Array of amounts [kmol] [kg] [kmol/kmol] [kg/kg]
    '''
    # Grabbing the total amount [kmol] and H2O content [kmol] for the mixture
    total = sum(mix.species_moles)
    gas = mix.phase_moles('gas')
    solid = mix.phase_moles('solid')
    H2Ocontent = mix.species_moles[pp.i['H2O']]
    # Pre-allocating the array
    amounts = np.zeros(len(species))
    # Getting the mole amounts [kmol] of the species
    for i, v in enumerate(species):
        amounts[i] = mix.species_moles[pp.i[v]]
    # Converting mole to mass [kg]
    if mass:
        for i, v in enumerate(species):
            amounts[i] *= pp.Mw[v]
        total = fs.getFuelMass(mix)
        # TODO: implement gas mass and solid mass calculations
        H2Ocontent *= pp.Mw['H2O']
    # Subtracting H2O content for dry basis
    if db:
        total -= H2Ocontent
        gas -= H2Ocontent
    # Normalizing the amounts
    if norm:
        if phase=='gas':
            amounts /= gas
        elif phase=='solid':
            amounts /= solid
        else:
            amounts /= total        
    # Making small values equal to zero
    if eps is not None:
        for j, a in enumerate(amounts):
            if a < eps:
                amounts[j] = 0
    return amounts

def H2CO(mix, mass=False):
    '''
    Returns the ratio of H2 moles to CO moles of the mixture.

    Parameters
    ----------
    mix : Cantera 'Mixture' object
        Fuel mixture object.
    mass : boolean
        If True, returns mass ratio [kg/kg]. (default: False, mole ratio [kmol/kmol])
    
    Returns
    -------
    H2CO : float
        Ratio of H2 amount to CO amount [kmol/kmol] [kg/kg]
    '''
    H2 = mix.species_moles[pp.i['H2']]
    CO = mix.species_moles[pp.i['CO']]
    if mass:
        H2 *= pp.Mw['H2']
        CO *= pp.Mw['CO']
    return H2/CO

def carbonConv(products, reagents):
    '''
    Returns the carbon conversion of the products to the reagents.

    Parameters
    ----------
    products : Cantera 'Mixture' object
        Object representing products.
    reagents : Cantera 'Mixture' object
        Object representing reagents.

    Returns
    -------
    conv : float
        Carbon conversion [dimensionless]
    '''
    Creag = sum(getAmounts(reagents, ['C(gr)', 'C']))
    Cprod = sum(getAmounts(products, ['C(gr)', 'C']))
    conv = (Creag - Cprod)/Creag
    return conv

def syngasHHV(mix, basis='vol', fuelMass = 1.0):
    '''
    Higher heating value of gas-phase products (syngas).

    Parameters
    ----------
    mix : Cantera 'Mixture' object
        Syngas mixture object.
    basis : string
        'mole' = mole amount in kmol
        'mass' = mass amount in kg
        'vol' = normal volume amount, Nm³
        Normal condition at 273.15 K and 1 atm.
        'fuel mass' = fuel mass in kg
    fuelMass : float
        Mass of fuel [kg] (default: 1.0)

    Returns
    -------
    HHV : float
        Higher heating value [MJ/kg] [MJ/kmol] [MJ/Nm³]
    '''
    HHV_H2 = pp.Hfo['H2'] + 0.5*pp.Hfo['O2'] - pp.Hfo['H2O']
    HHV_CH4 = pp.Hfo['CH4'] + 2*pp.Hfo['O2'] - pp.Hfo['CO2'] - 2*pp.Hfo['H2O']
    HHV_CO = pp.Hfo['CO'] + 1*pp.Hfo['O2'] - pp.Hfo['CO2']       # TODO: CHECK!!! 1 O2??? Shouldn't it be 0.5 O2?
    #HHV_C2H6 = pp.Hfo['C2H6'] + 3.5*pp.Hfo['O2'] - 2*pp.Hfo['CO2'] - 3*pp.Hfo['H2O']

    H2 = mix.species_moles[pp.i['H2']]
    CH4 = mix.species_moles[pp.i['CH4']]
    CO = mix.species_moles[pp.i['CO']]
    #C2H6 = mix.species_moles[pp.i['C2H6']]

    HHV = (1e-6)*(H2*HHV_H2 + CH4*HHV_CH4 + CO*HHV_CO) # + C2H6*HHV_C2H6)
    
    if basis == 'fuel mass':
        HHV /= fuelMass
    else:
        HHV /= gasYield(mix, basis)

    return HHV

def syngasLHV(mix):
    '''
    Lower heating value of gas-phase products (syngas).

    Parameters
    ----------
    mix : Cantera 'Mixture' object
        Syngas mixture object.

    Returns
    -------
    LHV : float
        Lower heating value [MJ/kmol]
    '''
    LHV_H2 = 120.092*pp.Mw['H2']
    LHV_CH4 = 49.855*pp.Mw['CH4']
    LHV_CO = 10.160*pp.Mw['CO']
    #LHV_C2H6 = 47.208*pp.Mw['C2H6']

    H2 = mix.species_moles[pp.i['H2']]
    CH4 = mix.species_moles[pp.i['CH4']]
    CO = mix.species_moles[pp.i['CO']]
    #C2H6 = mix.species_moles[pp.i['C2H6']

    LHV = H2*LHV_H2 + CH4*LHV_CH4 + CO*LHV_CO #+ C2H6*LHV_C2H6

    LHV *= (1 - mix.species_moles[pp.i['H2O']] / gasYield(mix, basis='mole', db=False)) # Não entendi essa linha de código

    # TODO: Add basis

    return LHV

# def gas_hhv(self, basis='vol'):
#     """
#     Higher heating value of gas-phase products (fuel gas).

#     Parameters
#     ----------
#     self : ndarray
#         Mole of products [kmol]
#     basis : string
#         HHV in mass fraction = 'w', mole fraction = 'm', 
#         volume fraction = 'v' (default)

#     Returns
#     -------
#     HHV : float
#         Higher heating value in the respective basis (mass, mole, or volume), 
#         d.b. [MJ/kg] [MJ/kmol] [MJ/Nm3]
#     """
#     ns = pp.nsp
#     # preallocate variables
#     sp = []
#     hhv_i = np.zeros(ns) # will be nonzero to 'heating' species
#     # find 'heating' species
#     for i in range(ns):
#         if (i == pp.f.species_index('gas','H2') or \
#             i == pp.f.species_index('gas','CO')):
#             # Combustion of hydrogen
#             # H2 + 0.5O2 --> H2O + <<HHV>>
#             # Combustion of carbon monoxide
#             # CO + 0.5O2 --> CO2 + <<HHV>>
#             sp = np.append(sp, pp.f.species_name(i))
#             hhv_i[i] = pp.Hfo[i] + (pp.f.n_atoms(i,'C') \
#             + 0.25*pp.f.n_atoms(i,'H'))*pp.Hfo[pp.i_O2] \
#             - (pp.f.n_atoms(i,'C'))*pp.Hfo[pp.i_CO2] \
#             # FIXME: liquid or gas water?
#             - (0.5*pp.f.n_atoms(i,'H'))*pp.Hfo[pp.i_H2O] # [J/kmol]
#         if (pp.f.n_atoms(i,'C') >= 1 and pp.f.n_atoms(i,'H') >= 1):
#             if (pp.f.n_atoms(i,'N') == 0 and pp.f.n_atoms(i,'O') == 0 and \
#                 pp.f.n_atoms(i,'S') == 0):
#                 # Combustion of hydrocarbons
#                 # CxHy + (x+0.25y)O2 --> xCO2 + 0.5yH2O + <<HHV>>
#                 sp = np.append(sp, pp.f.species_name(i))
#                 hhv_i[i] = pp.Hfo[i] + (pp.f.n_atoms(i,'C') \
#                 + 0.25*pp.f.n_atoms(i,'H'))*pp.Hfo[pp.i_O2] \
#                 - (pp.f.n_atoms(i,'C'))*pp.Hfo[pp.i_CO2] \
#                 # FIXME: liquid or gas water?
#                 - (0.5*pp.f.n_atoms(i,'H'))*pp.Hfo[pp.i_H2O] # [J/kmol]
#     ## N2 H2 CO CH4 CO2 C2H6
#     # higher heating value
#     hhv = np.sum(self*hhv_i)*1e-6 # [MJ]
#     if (basis == 'mole'):
#         return hhv/gas_yield(self, db='y', basis='mole') # d.b. [MJ/kmol]
#     if (basis == 'mass'):
#         return hhv/gas_yield(self, db='y', basis='mass') # d.b. [MJ/kg]
#     if (basis == 'vol'):
#         return hhv/gas_yield(self, db='y', basis='vol') # d.b. [MJ/Nm3]


def coldGasEff(mix, fuelLHV, moist=0.0):
    """
    Returns cold gas efficiency of gasification.

    Parameters
    ----------
    mix : Cantera 'Mixture' object
        Fuel mixture object.
    fuelLHV : float
        Fuel LHV
    moist : float
        Fuel moisture

    Returns
    -------
    coldGasEff : float
        Cold gas efficiency [dimensionless]
    """
    return syngasLHV(mix)/fuelLHV


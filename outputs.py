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
import pp as ppold
import feedstock as fsold
import pp2 as pp
import feedstock2 as fs
import cantera as ct
import numpy as np
import scipy.optimize as opt

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
    molSp = outlet.species_moles
    mass = fs.getFuelMass(outlet)

    if db:
        moles -= molSp[pp.i['H2O']]
        mass -= molSp[pp.i['H2O']] * pp.Mw['H2O']

    # wet basis       
    if basis == 'mole':
        return moles - molSp[pp.i['N2']]
    elif basis == 'mass':
        return mass - molSp[pp.i['N2']]*pp.Mw['N2']
    elif basis == 'vol':
        return ((moles - molSp[pp.i['N2']]))*R*Tn/Pn
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

def getAmounts(mix, species, norm=False, db=False, eps=None, mass=False):
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

    Returns
    -------
    amounts : array
        Array of amounts [kmol] [kg] [kmol/kmol] [kg/kg]
    '''
    # Grabbing the total amount [kmol] and H2O content [kmol] for the mixture
    total = sum(mix.species_moles)
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
        H2Ocontent *= pp.Mw['H2O']
    # Subtracting H2O content for dry basis
    if db:
        total -= H2Ocontent
    # Normalizing the amounts
    if norm:
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

# def syngas_hhv(self, fuel_mass=1.0, basis='vol'):
#     """
#     Higher heating value of gas-phase products (syngas).

#     Parameters
#     ----------
#     self : ndarray
#         Mole of products [kmol]
#     fuel : float
#         Mass of fuel, w.b.
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
#     # find key species
#     for i in range(ns):
#         if (i == pp.f.species_index('gas','H2') or \
#             i == pp.f.species_index('gas','CH4') or \
#             i == pp.f.species_index('gas','CO') #or \
# #            i == pp.f.species_index('gas','C2H6')
#             ):
#             sp = np.append(sp, pp.f.species_name(i))
#             hhv_i[i] = pp.Hfo[i] + (pp.f.n_atoms(i,'C') \
#             + 0.25*pp.f.n_atoms(i,'H'))*pp.Hfo[pp.i_O2] \
#             - (pp.f.n_atoms(i,'C'))*pp.Hfo[pp.i_CO2] \
#             # FIXME: liquid or gas water?
#             - (0.5*pp.f.n_atoms(i,'H'))*pp.Hfo[pp.i_H2O] # [J/kmol]
#     # higher heating value
#     hhv = np.sum(self*hhv_i)*1e-6 # [MJ]
#     if (basis == 'syngas mole'):
#         return hhv/gas_yield(self, db='y', basis='mole') # d.b. [MJ/kmol]
#     if (basis == 'syngas mass'):
#         return hhv/gas_yield(self, db='y', basis='mass') # d.b. [MJ/kg]
#     if (basis == 'fuel mass'):
#         return hhv/fuel_mass # [MJ/kg]
#     if (basis == 'syngas vol'):
#         return hhv/gas_yield(self, db='y', basis='vol') # d.b. [MJ/Nm3]

# def syngas_lhv(self, fuel_mass=1.0):
#     """
#     Lower heating value (LHV) of gas-phase products (syngas).

#     Parameters
#     ----------
#     self : ndarray
#         Mole of products [kmol]
#     fuel : float
#         Mass of fuel, w.b.
#     basis : string
#         LHV in mass fraction = 'w', mole fraction = 'm', 
#         volume fraction = 'v' (default)

#     Returns
#     -------
#     lhv : float
#         Lower heating value [MJ/kg]
#     """
#     lhv_CO = 10.160*pp.Mw[pp.i_CO] # MJ/kmol
#     lhv_CH4 = 49.855*pp.Mw[pp.i_CH4] # MJ/kmol
# #    lhv_C2H6 = 47.208*pp.Mw[pp.i_C2H6] # MJ/kmol
#     lhv_H2 = 120.092*pp.Mw[pp.i_H2] # MJ/kmol
#     return (lhv_CO*self[pp.i_CO] + lhv_CH4*self[pp.i_CH4] \
# #            + lhv_C2H6*self[pp.i_C2H6] 
#             + lhv_H2*self[pp.i_H2])*(1 \
#             - self[pp.i_H2O]/gas_yield(self, db='n', basis='mole'))
        
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

# def cold_gas_efficiency(self, fuel_lhv, moisture_fuel):
#     """
#     Return cold gas efficiency of gasification.
#     fuel_lhv : ndarray
#         Fuel LHV
#     moisture_fuel : ndarray
#         Fuel moisture
#     """
#     return (syngas_lhv(self, 1 + moisture_fuel)/fuel_lhv)    

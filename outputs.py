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

def gasYield(outlet, basis='vol', db='y'):
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

    if db == 'y':
        moles -= molSp[pp.i['H2O']]
        mass -= molSp[pp.i['H2O']] * pp.Mw['H2O']
    elif db != 'n':
        raise ValueError('db must be y or n')

    # wet basis       
    if basis == 'mole':
        return moles - molSp[pp.i['N2']]
    elif basis == 'mass':
        return mass - molSp[pp.i['N2']]*pp.Mw['N2']
    elif basis == 'vol':
        return ((moles - molSp[pp.i['N2']]))*R*Tn/Pn
    else:
        raise ValueError('Basis must be mole, mass or vol')

# def get_species(self, species=[], eps=1e-6):
#     '''
#     Get a list of species which mole fractions in 'self' are higher than 'eps'.
#     this function is useful to find a minimum number of species to handle out a
#     chemical equilibrium problem.
#     '''
#     i = 1
#     while i < pp.nsp:
#         if self[i] > eps:
#             species_name = pp.f.species_name(i)
#             try:
#                 species.index(species_name)
#             except:
#                 # exclude liquid species
#                 if 'L)' not in species_name:
#                     species.append(species_name)
#         i += 1
#     return species
    
# def get_fraction(self, species, normalized='n', db='n', eps=None):
#     '''
#     db : string
#         Dry basis ('y') or wet basis ('n', default)
#     '''
#     ## TODO: Make available for mass fraction calculation
#     idx = len(species)
#     mole = np.zeros(idx, 'd')
#     i = 0
#     while i < idx:
#         # get values
#         try:
#             mole[i] = self[pp.f.species_index('solid', species[i])]#/mole_solid
#         except:
#             mole[i] = self[pp.f.species_index('gas', species[i])]#/mole_gas
#         if eps != None:
#             # make small values as zero
#             if mole[i] < eps:
#                 mole[i] = 0
#         i += 1
#     # convert mole amount to mole fraction
#     mole /= sum(self)
#     if db == 'y':
#         mole *= (1 - self[pp.i_H2O])
#     if normalized == 'y':
#         # normalize values
#         mole /= np.sum(mole)
#     return mole
    
# def h2co_ratio(self):
#     h2 = self[pp.f.species_index('gas', 'H2')]
#     co = self[pp.f.species_index('gas', 'CO')]
#     return h2/co
    
# def carbon_conversion(products, reagents):
#     return (reagents[pp.i_C] - products[pp.i_C]) / reagents[pp.i_C]

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

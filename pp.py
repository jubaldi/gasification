#!/usr/bin/env python

"""This script sets thermophysical properties of phases and species.

The definition of phases must be done here. This file is necessary to 
other scripts.

@author = Rodolfo Rodrigues
@contact = rodolfo.rodrigues@ufsm.br
@data = April, 2012, rev.: June, 2013 (adapted to use cython Cantera)
"""
#==============================================================================
# import libraries
#==============================================================================
import sys
import cantera as ct
import numpy as np
import scipy.optimize as opt

import feedstock as fs
import fuel as fu
import gasifier as g
import outputs as op
import energy as en
#==============================================================================
# set parameters
#==============================================================================
To = 298.195 # K (reference temperature)
R = ct.gas_constant # 8314.47215 Pa*m^3/K/kmol

#==============================================================================
# load components to phases
#==============================================================================
s = ct.Solution('data/min-species.xml','solid')
g = ct.Solution('data/min-species.xml','gas')
#print s.species_names
air = ct.Solution('data/air.cti','air')
def mix():
    return ct.Mixture([(s,1),(g,1)])
f = ct.Mixture([(s,1),(g,1)])
# preallocate variable
nsp = f.n_species

airMix = ct.Mixture([(air,1)])

# get species names as lists
names_s = s.species_names
names_g = g.species_names
names_f = f.species_names

#==============================================================================
# get standard state properties
#==============================================================================
# get molecular weights as lists
Mw_s = s.molecular_weights
Mw_g = g.molecular_weights

# turn molecular weights into dict
Mw = {**dict(zip(names_s,Mw_s)), **dict(zip(names_g,Mw_g))}
#Mw['H2O(l)'] = Mw['H2O']
# air molecular weight (kg/kmol)
Mw_air = air.mean_molecular_weight

# standard enthalpy of formation
# get Hfo as lists
Hfo_s = s.standard_enthalpies_RT*To*R
Hfo_g = g.standard_enthalpies_RT*To*R
# turn Hfo into dict
Hfo = {**dict(zip(names_s,Hfo_s)), **dict(zip(names_g,Hfo_g))}
Hfo['H2O(l)'] = -283970115.359

## find important indices
# fuel species
# turn indices into dict
idx_dict_s = {name : f.species_index('solid',name) for name in names_s}
idx_dict_g = {name : f.species_index('gas',name) for name in names_g}
i = {**idx_dict_s, **idx_dict_g}
#i['H2O(l)'] = f.n_species
# create a new dict that inverts keys and values
i_inv = {v: k for k, v in i.items()}

# water enthalpy of vaporization
H_vap     = 4.0667e7 # J/kmol
#H_vap_mass= H_vap/Mw[i_H2Ol] # J/kg
H_vap_mass= H_vap/Mw['H2O'] # J/kg
Hfo_air   = 0.23211606*Hfo['O2'] + 0.75507754*Hfo['N2'] + 0.0128064*Hfo['Ar']

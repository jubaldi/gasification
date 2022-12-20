"""Instatiates phase objects and determines its properties.

This script sets thermophysical properties of phases, species and mixtures of phases.
All Cantera phase and mixture objects are instatiated here.

First made by: Rodolfo Rodrigues, rodolfo.rodrigues@ufsm.br
Revamp by: Nícolas Anése, nicolasanese19@gmail.com

Date: April, 2012, rev.: December, 2022
"""

# Imports libraries
import cantera as ct

# Sets parameters
To = 298.15 # K (reference temperature)
Po = 101325 # Pa (standard pressure)
R = ct.gas_constant # 8314.47215 Pa*m^3/K/kmol

# Instantiates phase (solution) objects with specific species
# according to xml file provided
solid = ct.Solution('data/min-species.xml','solid')
gas = ct.Solution('data/min-species.xml','gas')
air = ct.Solution('data/air.cti','air')

# Instantiates a mixture of phases with 0 moles of each species
def mix():
    m = ct.Mixture([(solid,0),(gas,0)])
    m.T = To
    m.P = Po
    return m

# New Fuel class inherits Mixture class from Cantera
# and adds some properties such as HHV, etc.
class Fuel(ct.Mixture):
    def testMethod(self):
        return [moles+1 for moles in self.species_moles]
    testAttribute = 10
    HHV = None
    LHV = None
    moisture = None
    # TO BE CONTINUED

def fuel():
    f = Fuel([(solid,0),(gas,0)])
    f.T = To
    f.P = Po
    return f

names_s = solid.species_names
names_g = gas.species_names

# get molecular weights as lists
Mw_s = solid.molecular_weights
Mw_g = gas.molecular_weights

# turn molecular weights into dict
Mw = {**dict(zip(names_s,Mw_s)), **dict(zip(names_g,Mw_g))}
Mw['H2O(l)'] = Mw['H2O']
    
def get_molecular_weight(species):
    return Mw[species]

# standard enthalpy of formation
# get Hfo as lists
Hfo_s = solid.standard_enthalpies_RT*To*R
Hfo_g = gas.standard_enthalpies_RT*To*R
# turn Hfo into dict
Hfo = {**dict(zip(names_s,Hfo_s)), **dict(zip(names_g,Hfo_g))}
Hfo['H2O(l)'] = -283970.115359 # J / kmol

def get_enthalpy_formation(species):
    return Hfo[species]
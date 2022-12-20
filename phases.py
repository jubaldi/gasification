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

# PROPERTIES - Molecular weight, Enthalpy of formation

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

# Important functions and classes

# Instantiates a mixture of phases with 0 moles of each species
def mix():
    m = ct.Mixture([(solid,0),(gas,0)])
    m.T = To
    m.P = Po
    return m

testMix = mix()
indices = {}
for index, species in enumerate(testMix.species_names):
    indices[species] = index

# New Fuel class inherits Mixture class from Cantera
# and adds some properties such as HHV, etc.
class Fuel(ct.Mixture):
    # def testMethod(self):
    #     return [moles+1 for moles in self.species_moles]
    # testAttribute = 10
    HHV = None
    LHV = None
    moisture = None
    ashFraction = None

    def get_mass(self):
        mass = 0
        for i, species in enumerate(self.species_names):
            speciesMole = self.species_moles[i]
            speciesMW = Mw[species]
            speciesMass = speciesMole*speciesMW
            mass += speciesMass
        return mass
    
    def set_moisture(self, moisture):
        # Sets moisture content to given fraction (%m/m dry basis)
        self.moisture = moisture
        mass = self.get_mass()
        currentWaterMoles = self.species_moles[self.species_index('gas', 'H2O')]
        currentWaterMass = currentWaterMoles * Mw['H2O']
        dryMass = mass - currentWaterMass
        newWaterMass = dryMass*moisture
        newWaterMoles = newWaterMass / Mw['H2O']
        newMoles = self.species_moles
        newMoles[indices['H2O']] = newWaterMoles
        self.species_moles = newMoles
    
    def redistribute_ash(self, ashComposition):
        # Sets a new ash composition, mantaining fixed ash mass, given ash composition distribution (%m/m of ash)
        ashComponents = ['SiO2(hqz)', 'CaO(s)', 'AL2O3(a)', 'Fe2O3(s)', 'Na2O(c)', 'K2O(s)', 'MgO(s)', 'P2O5', 'TiO2(ru)', 'SO3', 'Cr2O3(s)']
        ashMass = 0
        for index, species in enumerate(ashComponents):
            speciesMoles = self.species_moles[indices[species]]
            speciesMass = speciesMoles * Mw[species]
            ashMass += speciesMass
        newAshMassDistr = [ashMass*comp for comp in ashComposition]
        newMoles = self.species_moles # list to be updated
        for index, species in enumerate(ashComponents):
            newSpeciesMass = newAshMassDistr[index]
            newSpeciesMoles = newSpeciesMass / Mw[species]
            newMoles[indices[species]] = newSpeciesMoles
        self.species_moles = newMoles # updates list

    # TO BE CONTINUED

def fuel():
    f = Fuel([(solid,0),(gas,0)])
    f.T = To
    f.P = Po
    return f


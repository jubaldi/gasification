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
Mw_air = 0.21*Mw['O2'] + 0.78*Mw['N2'] + 0.01*Mw['Ar']
    
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

# New class "Stream" inherits Mixture class from Cantera
# and adds some properties such as HHV, etc.
# Also adds important methods to calculate gasification and log results
class Stream(ct.Mixture):
    '''
    This is a general-purpose class for representing fuel, gasifying agent and syngas streams.
    It is inherited from Cantera 'Mixture' class.
    The function 'stream()' is a quick way to instantiate an empty stream.
    Some attributes and methods are best applied to a specific type of stream (fuel, agent, outlet).
    When functions are applied to Stream objects, some attributes are used to register its history:
    for instance, gasifier.gasify_isot function registers fuelHHV on outlet because this information
    is important to determine cold gas efficiency.
    '''
    # Fuel-specific attributes
    fuelHHV = None
    fuelLHV = None
    fuelMoisture = None
    fuelAshFraction = None
    fuelDryMass = None

    # Outlet-specific attributes
    carbonConversion = None

    # Method for computing total mass of stream
    def get_mass(self):
        mass = 0
        for i, species in enumerate(self.species_names):
            speciesMole = self.species_moles[i]
            speciesMW = Mw[species]
            speciesMass = speciesMole*speciesMW
            mass += speciesMass
        return mass
    
    # Fuel-specific method for setting moisture content when creating a fuel stream (used in feedstock.py)
    def set_moisture(self, moisture):
        # Sets moisture content to given fraction (%m/m dry basis)
        self.fuelMoisture = moisture
        mass = self.get_mass()
        currentWaterMoles = self.species_moles[self.species_index('gas', 'H2O')]
        currentWaterMass = currentWaterMoles * Mw['H2O']
        dryMass = mass - currentWaterMass
        newWaterMass = dryMass*moisture
        newWaterMoles = newWaterMass / Mw['H2O']
        newMoles = self.species_moles
        newMoles[indices['H2O']] = newWaterMoles
        self.species_moles = newMoles
    
    # Fuel-specific method for setting ash composition when creating a fuel stream (used in feedstock.py)
    def redistribute_ash(self, ashComposition):
        # Sets a new ash composition, mantaining fixed ash mass, given ash composition distribution (%m/m of ash)
        ashComponents = ['SiO2(hqz)', 'CaO(s)', 'AL2O3(a)', 'Fe2O3(s)', 'Na2O(c)', 'K2O(s)', 'MgO(s)', 'P2O5', 'TiO2(ru)', 'SO3', 'Cr2O3(s)']
        ashMass = 0
        # Normalize ash composition
        totalAsh = sum(ashComposition)
        ashComposition = [ashComp / totalAsh for ashComp in ashComposition]

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
    
    def get_fuel_formula(self):
        nC = self.species_moles[indices['C(gr)']]
        nH = self.species_moles[indices['H']]
        nO = self.species_moles[indices['O']]
        nN = self.species_moles[indices['N']]
        nS = self.species_moles[indices['S']]
        nCl = self.species_moles[indices['CL']]
        nAsh = 0
        for i, comp in enumerate(['SiO2(hqz)', 'CaO(s)', 'AL2O3(a)', 'Fe2O3(s)', 'Na2O(c)', 'K2O(s)', 'MgO(s)', 'P2O5', 'TiO2(ru)', 'SO3', 'Cr2O3(s)']):
            nAsh += self.species_moles[indices[comp]]
        formula = [1, nH/nC, nO/nC, nN/nC, nS/nC, nCl/nC, nAsh/nC]
        return formula

    # Method for computing mole fraction of a certain species in the gas phase.
    def get_gas_fraction(self, species):
        # 'species' must be gas-phase species.
        gasMoles = self.phase_moles('gas')
        speciesMoles = self.species_moles[self.species_index('gas', species)]
        return speciesMoles / gasMoles

    # Method for computing amount of syngas. This can then be divided by fuelMass to obtain syngas yield.
    # Preferably used on outlet streams to determine syngas yield.
    def get_gas_amount(self, quantity='vol', basis='db'):
        
        gasAmount = None

        return gasAmount

    def get_syngas_hhv(self):
        HHV = None
        return HHV

# Instantiates a stream with 0 moles of each species
def stream():
    st = Stream([(solid,0),(gas,0)])
    st.T = To
    st.P = Po
    return st

# 'indices' is a dictionary that relates species names to their index.
testMix = stream()
indices = {}
for index, species in enumerate(testMix.species_names):
    indices[species] = index
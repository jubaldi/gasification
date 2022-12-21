"""This script defines functions to equilibrium simulation of gasification 
processes. It uses some predefined functions from Cantera package.

First made by: Rodolfo Rodrigues, rodolfo.rodrigues@ufsm.br
Revamp by: Nícolas Anése, nicolasanese19@gmail.com
Date: April, 2012, rev.: December, 2022
"""

#==============================================================================
# import libraries/files
#==============================================================================
import sys
import cantera as ct
import numpy as np
import scipy.optimize as opt
import matplotlib.pyplot as plt

import phases
import feedstock as fs
import fuel as fu
# import outputs as op
# import energy as en
#==============================================================================
# predefine parameters
#==============================================================================
R = ct.gas_constant  # 8314.4621 Pa*m^3/K/kmol
Tn = 273.15  # K
Pn = ct.one_atm  # 101315 Pa
zero = np.zeros(1)
one = np.ones(1)

def gasify_isot(fuel, agent, T=298.15, P=101325, carbonConv=1.0, directMethaneConv=0.0):
    # CARBON AVAILABILITY AND DIRECT METHANE CONVERSION
    # Procedure according to Li et al. (2004) (10.1016/S0961-9534(03)00084-9)
    # Fraction of carbon that goes through equilibrium process:
    carbonEquilibriumFraction = carbonConv - directMethaneConv
    carbonMoles = fuel.species_moles[phases.indices['C(gr)']]
    carbonEquilibriumMoles = carbonMoles * carbonEquilibriumFraction
    directMethaneMoles = carbonMoles * directMethaneConv
    charMoles = carbonMoles * (1 - carbonConv)
    consumedHydrogenMoles = 4 * directMethaneMoles
    newMoles = fuel.species_moles
    newMoles[phases.indices['C(gr)']] = carbonEquilibriumMoles
    newMoles[phases.indices['H']] = newMoles[phases.indices['H']] - consumedHydrogenMoles
    
    outlet = phases.stream()

    outletMoles = [fuelMole + agent.species_moles[i] for i, fuelMole in enumerate(newMoles)]
    outlet.species_moles = outletMoles
    outlet.T = T
    outlet.P = P
    outlet.equilibrate('TP')

    # At the end, char and methane are added back to the stream
    outletMoles = outlet.species_moles
    outletMoles[phases.indices['C(gr)']] += charMoles
    outletMoles[phases.indices['CH4']] += directMethaneMoles
    outlet.species_moles = outletMoles

    outlet.fuelMoisture = fuel.fuelMoisture
    outlet.fuelHHV = fuel.fuelHHV
    outlet.fuelLHV = fuel.fuelLHV
    outlet.fuelAshFraction = fuel.fuelAshFraction

    return outlet

def gasify(fuel, agent, T0=298.15, P=101325, carbonConv=1.0, directMethaneConv=0.0, isot=False):
    pass



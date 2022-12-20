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

def gasify_isot(fuel, agent, T, P):
    outlet = phases.outlet()

    outletMoles = [fuelMole + agent.species_moles[i] for i, fuelMole in enumerate(fuel.species_moles)]
    outlet.species_moles = outletMoles
    outlet.T = T
    outlet.P = P
    outlet.equilibrate('TP')

    outlet.fuelMoisture = fuel.fuelMoisture
    outlet.fuelHHV = fuel.fuelHHV
    outlet.fuelLHV = fuel.fuelLHV
    outlet.fuelAshFraction = fuel.fuelAshFraction

    return outlet

fuel1 = fs.create_fuel_stream(1, [57.2/91, 3.3/91, 16.2/91, 0.7/91, 0.2/91, 0], 13.4/91, 9/91, 21.1)
air1 = fs.create_air_stream(10, 298.15, 101325)
outlet1 = gasify_isot(fuel1, air1, 700, 101325)
print(fs.convert_ER_air(fuel1, 0.5))
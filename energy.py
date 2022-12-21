
"""
Defines functions to calculate energy parameters required for
non-isothermic equilibrium simulation of gasification 
processes. It uses some predefined functions from Cantera package.

First made by: Rodolfo Rodrigues, rodolfo.rodrigues@ufsm.br
Revamp by: Nícolas Anése, nicolasanese19@gmail.com
Date: April, 2012, rev.: December, 2022
"""

import sys
import cantera as ct
import numpy as np
import scipy.optimize as opt 

import phases
import feedstock as fs

R = ct.gas_constant  # 8314.4621 Pa*m^3/K/kmol
Tn = 273.15  # K
Pn = ct.one_atm  # 101315 Pa
zero = np.zeros(1)
one = np.ones(1)

def get_fuel_hf_Li(fuel):
    HHV = fuel.fuelHHV*1000 # kJ/kg
    moistureDB = fuel.fuelMoisture
    moistureWB = moistureDB / (1 + moistureDB)
    fuelMass = fuel.get_mass()
    C_ar = fuel.species_moles[phases.indices['C(gr)']] * phases.Mw['C(gr)'] / fuelMass
    H_ar = fuel.species_moles[phases.indices['H']] * phases.Mw['H'] / fuelMass
    S_ar = fuel.species_moles[phases.indices['S']] * phases.Mw['S'] / fuelMass
    # Correlation by Li et al. (2001) (10.1016/S0016-2361(00)00074-0)
    DHfo = HHV - (327.63*C_ar + 1417.94*H_ar + 92.57*S_ar + 158.67*moistureWB) # kJ/kg
    return DHfo # kJ/kg

ultimateDAF = [50.61, 8.79, 25.46, 12.80, 1.89, 0.46]
ashWB = 0.051
moistWB = 0.124
ashDB = ashWB / (1 - moistWB)
moistDB = moistWB / (1 - moistWB)
leather = fs.create_fuel_stream(1, ultimateDAF, ashDB, moistDB, 18.448, 16.729, 
                                [23.58, 1.20, 7.35, 2.59, 1.08, 0.79, 0, 1.59, 1.53, 0, 55.91])
emptyAir = fs.create_air_stream(fs.convert_ER_air(leather, 0))

print(get_fuel_hf_Li(leather))
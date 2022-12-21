
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

# def get_fuel_hf_Li(fuel):
#     HHV = fuel.fuelHHV*1000 # kJ/kg
#     moistureDB = fuel.fuelMoisture
#     moistureWB = moistureDB / (1 + moistureDB)
#     fuelMass = fuel.get_mass()
#     C_ar = fuel.species_moles[phases.indices['C(gr)']] * phases.Mw['C(gr)'] / fuelMass
#     H_ar = fuel.species_moles[phases.indices['H']] * phases.Mw['H'] / fuelMass
#     S_ar = fuel.species_moles[phases.indices['S']] * phases.Mw['S'] / fuelMass
#     # Correlation by Li et al. (2001) (10.1016/S0016-2361(00)00074-0)
#     DHfo = HHV - (327.63*C_ar + 1417.94*H_ar + 92.57*S_ar + 158.67*moistureWB) # kJ/kg
#     return DHfo # kJ/kg

def get_fuel_enthalpy_formation(fuel):
    '''
    Estimates the standard enthalpy of formation of the given fuel [J]
    from Higher Heating Value (HHV) and species composition.

    Parameters
    ----------
    fuel : Fuel object
        Fuel object

    Returns
    -------
    hFormation : float
        Standard enthalpy of formation [J]
    '''
    h = lambda x: phases.Hfo[x]
    n = lambda x: fuel.species_moles[phases.indices[x]]
    HHV_by_fuelMass = fuel.fuelHHV # MJ/kg
    HHV = HHV_by_fuelMass * fuel.fuelDryMass * 1E6 # J

    stoic = fs.determine_stoich(fuel) * sum(fuel.species_moles)
    hFormation1 = (n('C(gr)')*h('CO2') + 0.5*n('H')*h('H2O(l)') + n('S')*h('SO2') + n('CL')*h('CLO') + HHV) / n('C(gr)')
    hFormationN = (n('N')*h('N2')) / n('C(gr)')
    hFormationO = (-stoic*h('O2')) / n('C(gr)')
    hFormationAsh = (n('CaO(s)')*h('CaO(s)') +                  
                    n('SiO2(hqz)')*h('SiO2(hqz)') + n('AL2O3(a)')*h('AL2O3(a)') +
                    n('Fe2O3(s)')*h('Fe2O3(s)') + n('Na2O(c)')*h('Na2O(c)') + 
                    n('K2O(s)')*h('K2O(s)') + n('MgO(s)')*h('MgO(s)') + 
                    n('P2O5')*h('P2O5') + n('TiO2(ru)')*h('TiO2(ru)') + 
                    n('SO3')*h('SO3') + n('Cr2O3(s)')*h('Cr2O3(s)')) / n('C(gr)')
    hFormation = hFormation1 + hFormationN + hFormationO + hFormationAsh # J / (1 kmol C(gr))
    hFormation_total = hFormation * fuel.species_moles[phases.indices['C(gr)']]
    hFormation_mass = hFormation_total / fuel.get_mass() # J/kg
    # EDITED 31-10-22: Removed division by moles of C, which I thought was dimensionally inconsistent (Nícolas)
    # EDITED 20-12-22: Went back on that division, turns out it was correct. Oops (Nícolas)
    return hFormation_total

def get_enthalpy(stream):
    # Returns stream enthalpy in J.
    # Reference state is T = 298.15 K.
    enthalpy = 0
    for i, phase in enumerate(['solid', 'gas']):
        if stream.phase_moles(i) == 0:
            enthalpy += 0
        else:
            enthalpy += stream.phase(i).enthalpy_mole * stream.phase_moles(i)
    return enthalpy

ultimateDAF = [50.61, 8.79, 25.46, 12.80, 1.89, 0.46]
ashWB = 0.051
moistWB = 0.124
ashDB = ashWB / (1 - moistWB)
moistDB = moistWB / (1 - moistWB)
leather = fs.create_fuel_stream(1, ultimateDAF, ashDB, moistDB, 18.448, 16.729, 
                                [23.58, 1.20, 7.35, 2.59, 1.08, 0.79, 0, 1.59, 1.53, 0, 55.91])
emptyAir = fs.create_air_from_ER(leather, 0)


formula = leather.get_fuel_formula()
air1 = fs.create_air_from_ER(leather, 0.5)
air1.T = 298.15
print(air1.phase(1).enthalpy_mole)
steam = fs.create_steam_stream(10)
print(air1.phase(1).enthalpy_mole)

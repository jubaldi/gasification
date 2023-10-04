
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
import matplotlib.pyplot as plt

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
    DHfo_total = DHfo * fuel.fuelDryMass * 1000 # J
    return DHfo_total # J

def get_dry_fuel_HF(fuel):
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
        Standard enthalpy of formation [J/kmol]
    '''
    h = lambda x: phases.Hfo[x]
    n = lambda x: fuel.species_moles[phases.indices[x]]
    HHV_by_fuelMass = fuel.fuelHHV # MJ/kg, dry basis
    HHV = HHV_by_fuelMass * fuel.get_dry_mass() * 1E6 # J

    stoic = fs.determine_stoich(fuel) * sum(fuel.species_moles)
    hFormationCombs = (n('C(gr)')*h('CO2') + 0.5*n('H')*h('H2O(l)') + n('S')*h('SO2') + n('CL')*h('CLO')) / n('C(gr)')
    hFormationHHV = HHV / n('C(gr)')
    # hFormationMoist = (n('H2O')*h('H2O(l)')) / n('C(gr)')
    # hFormationN = (n('N')*h('N2')) / n('C(gr)')
    # hFormationO = (-stoic*h('O2')) / n('C(gr)')
    ashMoles = fuel.fuelDryMass*fuel.fuelAshFraction/fuel.fuelAshMW
    hFormationAsh = fuel.fuelAshHF * ashMoles / n('C(gr)')
    # print('Combs', hFormationCombs*1E-6)
    # print('HHV', hFormationHHV*1E-6)
    # print('Ash', hFormationAsh*1E-6)
    hFormation = hFormationCombs + hFormationHHV + hFormationAsh # J / (1 kmol C(gr))
    hFormation_total = hFormation * fuel.species_moles[phases.indices['C(gr)']]
    hFormation_mass = hFormation_total / fuel.get_mass() # J/kg
    return hFormation

# coalUltimate = [78.4750, 3.9681, 16.0249, 0.7044, 0.7748, 0.0528] # % daf
# coalAshDB = 0.130816 # fraction, w.b.
# coalMoistDB = 0 # fraction, w.b.
# coalHHV = 25.232 # MJ/kg
# coalLHV = 24.648 # MJ/kg
# ashComposition = [54.06, 6.57, 23.18, 6.85, 0.82, 1.6, 1.83, 0.5, 1.05, 3.54, 0] # % of ash
# coal = fs.create_fuel_stream(980, coalUltimate, coalAshDB, coalMoistDB, HHV=coalHHV, LHV=coalLHV, ashComposition=ashComposition)
# coal.fuelAshMW = 80.7
# coal.fuelAshHF = -788.92*1E6
# print(get_dry_fuel_HF(coal)*1E-6)

def get_enthalpy(stream):
    # Returns stream enthalpy in J.
    # Reference state is T = 298.15 K.
    enthalpy = 0
    for i, phase in enumerate(['solid', 'gas']):
        stream.phase(i).basis = 'mass'
        if stream.phase_moles(i) == 0:
            enthalpy += 0
        else:
            enthalpy += stream.phase(i).h * stream.get_phase_mass(i)
    return enthalpy

def get_enthalpy_formation(stream):
    # Returns enthalpy of formation for stream.
    enthalpyFormation = 0
    for i, sp in enumerate(stream.species_names):
        enthalpyFormation += phases.Hfo[sp] * stream.species_moles[i]
    return enthalpyFormation



# ultimateDAF = [50.61, 8.79, 25.46, 12.80, 1.89, 0.46]
# ashWB = 0.051
# moistWB = 0.124
# ashDB = ashWB / (1 - moistWB)
# moistDB = moistWB / (1 - moistWB)
# leather = fs.create_fuel_stream(1, ultimateDAF, ashDB, moistDB, 18.448, 16.729, 
#                                 [23.58, 1.20, 7.35, 2.59, 1.08, 0.79, 0, 1.59, 1.53, 0, 55.91])

# methane = fs.create_fuel_stream(1, [74.8705303, 25.1294697, 0, 0, 0, 0], 0, 0, 55.51954123)
# moles = methane.species_moles
# moles[phases.indices['C(gr)']] = 1
# moles[phases.indices['H']] = 4
# methane.species_moles = moles
# print(methane.species_moles[phases.indices['C(gr)']])
# print(methane.species_moles[phases.indices['H']])
# h1 = get_fuel_enthalpy_formation(methane)
# print(fs.determine_stoich(methane))
# air = fs.create_air_from_ER(leather, 0.2)
# h2 = get_enthalpy(air)
# h3 = get_enthalpy_formation(air)
# h0 = h1 + h2 + h3
# outlet = phases.stream()
# print(len(outlet.species_names))

# outletMoles = [fuelMole + air.species_moles[i] for i, fuelMole in enumerate(leather.species_moles)]
# outlet.species_moles = outletMoles
# outlet.T = 298.15
# outlet.P = 101325
# outlet.equilibrate('TP')
# TT = np.arange(200, 3200+1E-6, 50)
# HHf = []
# HHp = []
# HH = []
# for i, t in enumerate(TT):
#     outlet.T = t
#     outlet.P = 101325
#     outlet.equilibrate('TP')
#     hf = get_enthalpy_formation(outlet)
#     hp = get_enthalpy(outlet)
#     HHf.append(hf)
#     HHp.append(hp)
#     HH.append(hf + hp)
    

# plt.plot(TT, HHp)
# plt.plot(TT, HHf)
# plt.plot(TT, HH)
# plt.axhline(y=h0)
# plt.axhline(y=0)
# plt.grid()
# plt.xlim(TT[0], TT[-1])
# plt.show()


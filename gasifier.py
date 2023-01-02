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
import energy as en
# import outputs as op
#==============================================================================
# predefine parameters
#==============================================================================
R = ct.gas_constant  # 8314.4621 Pa*m^3/K/kmol
Tn = 273.15  # K
Pn = ct.one_atm  # 101315 Pa
zero = np.zeros(1)
one = np.ones(1)

def gasify_isot(fuel, agent, T=298.15, P=101325, charFormation=0.0, directMethaneConv=0.0):
    # CHAR FORMATION AND DIRECT METHANE CONVERSION
    # Procedure according to Li et al. (2004) (10.1016/S0961-9534(03)00084-9)
    # Fraction of carbon that goes through equilibrium process:
    carbonEquilibriumFraction = 1 - charFormation - directMethaneConv
    carbonMoles = fuel.species_moles[phases.indices['C(gr)']]
    carbonEquilibriumMoles = carbonMoles * carbonEquilibriumFraction
    directMethaneMoles = carbonMoles * directMethaneConv
    charMoles = carbonMoles * charFormation
    consumedHydrogenMoles = 4 * directMethaneMoles
    newMoles = fuel.species_moles
    newMoles[phases.indices['C(gr)']] = carbonEquilibriumMoles
    newMoles[phases.indices['H']] = newMoles[phases.indices['H']] - consumedHydrogenMoles
    
    outlet = phases.stream()

    outletMoles = [fuelMole + agent.species_moles[i] for i, fuelMole in enumerate(newMoles)]
    outlet.species_moles = outletMoles
    outlet.T = T
    outlet.P = P
    outlet.equilibrate('TP', solver='gibbs', max_steps=1E6)

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

def gasify_nonisot(fuel, agent, T0=298.15, P=101325, heatLossFraction=0.0, charFormation=0.0, directMethaneConv=0.0, isot=False):
    fuel.T = T0
    fuel.P = P
    agent.T = T0
    agent.P = P
    enthalpyFormationFuel = en.get_dry_fuel_HF(fuel) * fuel.species_moles[phases.indices['C(gr)']] # J
    enthalpyFormationMoisture = (fuel.fuelDryMass * fuel.fuelMoisture / phases.Mw['H2O(l)']) * phases.Hfo['H2O(l)'] # J
    enthalpyAgent = en.get_enthalpy(agent) # J
    # enthalpyFormationAgent = en.get_enthalpy_formation(agent)
    initialEnthalpy = enthalpyFormationFuel + enthalpyFormationMoisture + enthalpyAgent # J
    lostEnthalpy = abs(enthalpyFormationFuel) * heatLossFraction
    desiredEnthalpy = initialEnthalpy - lostEnthalpy
    # print('formation fuel: ', enthalpyFormationFuel)
    # print('agent: ', enthalpyAgent)
    # # print('formation agent: ', enthalpyFormationAgent)
    # print('initial: ', initialEnthalpy)

    # CHAR FORMATION AND DIRECT METHANE CONVERSION
    # Procedure according to Li et al. (2004) (10.1016/S0961-9534(03)00084-9)
    # Fraction of carbon that goes through equilibrium process:
    carbonEquilibriumFraction = 1 - charFormation - directMethaneConv
    carbonMoles = fuel.species_moles[phases.indices['C(gr)']]
    carbonEquilibriumMoles = carbonMoles * carbonEquilibriumFraction
    directMethaneMoles = carbonMoles * directMethaneConv
    charMoles = carbonMoles * charFormation
    consumedHydrogenMoles = 4 * directMethaneMoles
    newMoles = fuel.species_moles
    newMoles[phases.indices['C(gr)']] = carbonEquilibriumMoles
    newMoles[phases.indices['H']] = newMoles[phases.indices['H']] - consumedHydrogenMoles

    outlet = phases.stream()

    outletMoles = [fuelMole + agent.species_moles[i] for i, fuelMole in enumerate(newMoles)]
    outlet.species_moles = outletMoles
    outlet.T = T0
    outlet.P = P
    outlet.equilibrate('TP')

    initialMoles = outlet.species_moles

    def residual(Teq):
        outlet.T = Teq
        outlet.P = P
        outlet.equilibrate('TP')
        # At the end, char and methane are added back to the stream
        outletMoles = initialMoles.copy()
        outletMoles[phases.indices['C(gr)']] += charMoles
        outletMoles[phases.indices['CH4']] += directMethaneMoles
        outlet.species_moles = outletMoles
        finalEnthalpy = en.get_enthalpy(outlet)
        error = (finalEnthalpy - desiredEnthalpy)**2
        return error
    
    solution = opt.minimize_scalar(residual, method='bounded', bounds=(200, 10000))
    Teq = solution.x

    outlet.species_moles = initialMoles
    outlet.T = Teq
    outlet.P = P
    outlet.equilibrate('TP')
    # At the end, char and methane are added back to the stream
    outletMoles = outlet.species_moles
    outletMoles[phases.indices['C(gr)']] += charMoles
    outletMoles[phases.indices['CH4']] += directMethaneMoles
    outlet.species_moles = outletMoles
    # print('last: ', en.get_enthalpy(outlet))

    return outlet

# def gasify_nonisot2(fuel, agent, T0=298.15, P=101325, heatLossFraction=0.0, charFormation=0.0, directMethaneConv=0.0, isot=False):
#     fuel.T = T0
#     fuel.P = P
#     agent.T = T0
#     agent.P = P
#     enthalpyFormationFuel = en.get_dry_fuel_HF(fuel) * fuel.species_moles[phases.indices['C(gr)']] # J
#     enthalpyFormationMoisture = (fuel.fuelDryMass * fuel.fuelMoisture / phases.Mw['H2O(l)']) * phases.Hfo['H2O(l)'] # J
#     enthalpyAgent = en.get_enthalpy(agent) # J
#     # enthalpyFormationAgent = en.get_enthalpy_formation(agent)
#     initialEnthalpy = enthalpyFormationFuel + enthalpyFormationMoisture + enthalpyAgent # J
#     lostEnthalpy = abs(enthalpyFormationFuel) * heatLossFraction
#     desiredEnthalpy = initialEnthalpy - lostEnthalpy
#     # print('formation fuel: ', enthalpyFormationFuel)
#     # print('agent: ', enthalpyAgent)
#     # # print('formation agent: ', enthalpyFormationAgent)
#     # print('initial: ', initialEnthalpy)
    
#     def residual(Teq):
#         outlet = gasify_isot(fuel, agent, T=Teq, P=P, charFormation=charFormation, directMethaneConv=directMethaneConv)
#         finalEnthalpy = en.get_enthalpy(outlet)
#         error = (finalEnthalpy - desiredEnthalpy)**2
#         return error
    
#     solution = opt.minimize_scalar(residual, method='bounded', bounds=(200, 10000))
#     Teq = solution.x

#     outlet = gasify_isot(fuel, agent, T=Teq, P=P, charFormation=charFormation, directMethaneConv=directMethaneConv)

#     return outlet

# coalUltimate = [78.4750, 3.9681, 16.0249, 0.7044, 0.7748, 0.0528] # % daf
# coalAshDB = 0.130816 # fraction, d.b.
# coalMoistDB = 0.02 / (1 - 0.02) # fraction, d.b.
# coalHHV = 25.232 # MJ/kg
# coalLHV = 24.648 # MJ/kg
# ashComposition = [54.06, 6.57, 23.18, 6.85, 0.82, 1.6, 1.83, 0.5, 1.05, 3.54, 0] # % of ash
# coal = fs.create_fuel_stream(980, coalUltimate, coalAshDB, coalMoistDB, HHV=coalHHV, LHV=coalLHV, ashComposition=ashComposition)
# coal.fuelAshMW = 80.7 # kg/kmol
# coal.fuelAshHF = -788.92*1E6 # MJ/kmol
# air = fs.create_air_from_ER(coal, 1)
# outlet = gasify_isot(coal, air, T=1273.15, charFormation=0.8676)
# print(outlet.get_phase_mass(0))
# print(outlet.get_phase_mass(1))


# oltenia = fs.create_fuel_stream(1, [57.67, 4.59, 19.95, 1.56, 2.83, 0], 0, 13.39/(100-13.39))
# oxygen = fs.create_O2_from_ER(oltenia, 0.35)
# # print(air2.get_mass())
# Teq = 580 + 273.15 # K
# P = 1E6
# outlet = gasify_isot(oltenia, oxygen, Teq, P)
# print(outlet.species_moles)
# print('CO = ', outlet.get_gas_fraction('CO')*100)
# print('CO2 = ',outlet.get_gas_fraction('CO2')*100)
# print('CH4 = ',outlet.get_gas_fraction('CH4')*100)
# print('H2 = ',outlet.get_gas_fraction('H2')*100)
# print('N2 = ',outlet.get_gas_fraction('N2')*100)
# print('H2O = ',outlet.get_gas_fraction('H2O')*100)

# centralia = fs.create_fuel_stream(1, [76, 4.7, 18.7, 0.1, 0.4, 0], 0, 0)
# air = fs.create_air_from_ER(centralia, 0.3)
# # print(air2.get_mass())
# Teq = 585 + 273.15 # K
# P = 0.125E6
# outlet2 = gasify_isot(centralia, air, Teq, P)
# print(outlet2.species_moles)
# print('CO = ', outlet2.get_gas_fraction('CO')*100)
# print('CO2 = ',outlet2.get_gas_fraction('CO2')*100)
# print('CH4 = ',outlet2.get_gas_fraction('CH4')*100)
# print('H2 = ',outlet2.get_gas_fraction('H2')*100)
# print('N2 = ',outlet2.get_gas_fraction('N2')*100)
# print('H2O = ',outlet2.get_gas_fraction('H2O')*100)
# print(outlet2.get_syngas_hhv())

# ultimateDAF = [50.61, 8.79, 25.46, 12.80, 1.89, 0.46]
# ashWB = 0.051
# moistWB = 0.124
# ashDB = ashWB / (1 - moistWB)
# moistDB = moistWB / (1 - moistWB)
# HHV = 18.448
# LHV = 16.729
# ashComp = [23.58, 1.20, 7.35, 2.59, 1.08, 0.79, 0, 1.59, 1.53, 0, 55.91]
# ultimateDAF = [59.30, 6.16, 29.63, 1.08, 3.82, 0.01]
# ashWB = 0.431
# moistWB = 0.111
# ashDB = ashWB / (1 - moistWB)
# moistDB = moistWB / (1 - moistWB)
# HHV = 13.386
# LHV = 12.686
# ashComp = [29.54, 1.02, 9.02, 4.21, 0.16, 0.70, 0.41, 0.03, 0.36, 0, 0]
# fuel = fs.create_fuel_stream(1, ultimateDAF, ashDB, 0, HHV, LHV, ashComp)
# ERs = np.arange(0.5, 1+1E-6, 1E-2)
# CO2 = []
# CO = []
# C = []
# H2O = []
# H2 = []
# CH4 = []
# for i, er in enumerate(ERs):
#     air = fs.create_air_from_ER(fuel, er)
#     outlet = gasify_nonisot(fuel, air, T0=300, P=1E5, heatLossFraction=0.0, charFormation=0.0, directMethaneConv=0.0)
#     co2 = outlet.get_gas_fraction('CO2')
#     co = outlet.get_gas_fraction('CO')
#     c = outlet.species_moles[phases.indices['C(gr)']] / sum(outlet.species_moles)
#     h2o = outlet.get_gas_fraction('H2O')
#     h2 = outlet.get_gas_fraction('H2')
#     ch4 = outlet.get_gas_fraction('CH4')
#     H2O.append(h2o*100)
#     CO2.append(co2*100)
#     CO.append(co*100)
#     C.append(c*100)
#     H2.append(h2*100)
#     CH4.append(ch4*100)


# plt.plot(ERs, CO2, label='CO2')
# plt.plot(ERs, CO, label='CO')
# plt.plot(ERs, C, label='C(s)')
# plt.plot(ERs, H2O, label='H2O')
# plt.plot(ERs, H2, label='H2')
# plt.plot(ERs, CH4, label='CH4')
# plt.legend()
# plt.xlim(0.5, 1)
# # plt.ylim(0, 1800)
# plt.grid()
# plt.show()

# specs = [] 
# mols = []
# totalMoles = sum(outlet.species_moles)
# for i, mole in enumerate(outlet.species_moles):
#     frac = mole / totalMoles
#     if frac > 1E-4:
#         mols.append(mole)
#         specs.append(outlet.species_names[i])
# print(len(specs))
# print(specs)
# for i, sp in enumerate(specs):
#     print(f'{sp} = {outlet.get_gas_fraction(sp)*100}')


# ultimateDAF = [50.61, 8.79, 25.46, 12.80, 1.89, 0.46]
# ashWB = 0.051
# moistWB = 0.124
# ashDB = ashWB / (1 - moistWB)
# moistDB = moistWB / (1 - moistWB)
# leather = fs.create_fuel_stream(1, ultimateDAF, ashDB, moistDB, 18.448, 16.729, 
#                                 [23.58, 1.20, 7.35, 2.59, 1.08, 0.79, 0, 1.59, 1.53, 0, 55.91])
# air = fs.create_air_from_ER(leather, 0.4)

# outlet = gasify_nonisot(leather, air, T0=700)
# print('Teq: ', outlet.T)
# print('Heq: ', en.get_enthalpy(outlet) + en.get_enthalpy_formation(outlet))
# print(max(outlet.species_moles))
# print(sum(outlet.species_moles))

# specs = [] 
# mols = []
# totalMoles = sum(outlet.species_moles)
# for i, mole in enumerate(outlet.species_moles):
#     frac = mole / totalMoles
#     if frac > 1E-9:
#         mols.append(mole)
#         specs.append(outlet.species_names[i])

# print(len(specs))
# print(specs)
# print(phases.testMix.species_names)

# test = phases.stream()
# test.species_moles = [1.34287155e-09, 4.38909030e-05, 1.30290816e-05, 9.87523712e-06,
#  5.10638793e-06, 0.00000000e+00, 1.06096182e-05, 2.38951368e-04,
#  1.16641277e-05, 2.23972189e-04, 2.33041253e-11, 6.19729501e-13,
#  3.61650072e-02, 3.51389679e-03, 1.07061029e-02, 2.45279767e-02,
#  2.57233731e-03, 1.76728437e-02, 6.87565234e-04, 1.96021316e-04,
#  1.77095202e-06, 8.32993951e-02, 4.14083921e-08, 2.19804813e-08,
#  4.77210897e-04, 1.05544942e-04, 9.30256936e-07, 4.72104475e-05,
#  5.32436036e-06, 2.60598044e-04, 1.34266005e-04, 3.91085354e-09,
#  3.34549469e-07, 7.35640047e-12, 2.79922510e-05, 4.63049525e-11,
#  9.41928847e-05, 8.24428631e-09, 4.15504189e-15, 7.53492967e-12,
#  2.62632352e-12, 8.80426340e-06, 4.83627031e-06, 5.11456650e-15,
#  7.35477009e-22, 4.41161221e-34, 3.48801167e-21, 1.01584892e-03]
# test.T = 3005.219018769902
# # print(test.species_moles)
# # print(test.T)
# print(en.get_enthalpy(test) + en.get_enthalpy_formation(test))


"""Defines functions to create inputs (feedstock) for the gasification simulations.

First made by: Rodolfo Rodrigues, rodolfo.rodrigues@ufsm.br
Revamp by: Nícolas Anése, nicolasanese19@gmail.com
Date: April, 2012, rev.: December, 2022
"""

import cantera as ct
import numpy as np
import scipy.optimize as opt

import phases
import fuel as fu

def create_fuel_from_db(fuelID, fuelDryMass, moisture):
    '''
    Create a Cantera 'Mixture' object representing the given fuel.

    Parameters
    ----------
    fuelID : str
        The fuel ID as given by the CSV.
    fuelDryMass : float
        The fuel mass [kg] in dry basis
    moisture : float
        The fuel moisture content (%m/m, dry basis)

    Returns
    -------
    fuelMix : Cantera 'Mixture' object
        The fuel mixture object.
    '''
    comp = fu.fuel_comp(fuelID)
    massBySpecies = {key: value*fuelDryMass for key, value in comp.items()} # kg

    molesBySpecies = {key: value/phases.get_molecular_weight(key) for key, value in massBySpecies.items()} # kmol

    fuelMix = phases.stream()
    fuelMix.fuelHHV = fu.get_hv(fuelID, 'HHV')
    fuelMix.fuelLHV = fu.get_hv(fuelID, 'LHV')

    molesInOrder = np.zeros(len(fuelMix.species_names))
    for index, species in enumerate(fuelMix.species_names):
        if species in molesBySpecies.keys():
            molesInOrder[index] = molesBySpecies[species]
        else:
            molesInOrder[index] = 0

    fuelMix.species_moles = molesInOrder
    fuelMix.fuelDryMass = fuelDryMass

    fuelMix.set_moisture(moisture)

    return fuelMix

def create_fuel_stream(fuelDryMass, ultimate, ashFraction, moisture, HHV=None, LHV=None, 
                ashComposition=[50, 20, 10, 10, 10, 0, 0, 0, 0, 0, 0]):
    '''
    Create a Cantera 'Mixture' object representing the given fuel.

    Parameters
    ----------
    fuelDryMass : float
        The fuel mass [kg] in dry basis
    ultimate : list
        Ultimate analysis of the fuel: [C H O N S Cl] (mass fraction dry ash free (DAF) basis)
    ashFraction : float
        Fraction of ash in the fuel (mass fraction dry basis)
    moisture : float
        Moisture content of the fuel (mass fraction dry basis)
    HHV : float
        Higher heating value of the fuel (MJ/kg)
    LHV : float
        Lower heating value of the fuel (MJ/kg)

    Returns
    -------
    fuelMix : Fuel object
        The fuel mixture object.
    '''
    fuelMix = phases.stream()

    # Normalization: ultimate should add up to 100% (dry ash free basis)
    total = sum(ultimate)
    ultimate = [item / total for item in ultimate]
    # Ash compositon should be normalized as well
    totalAsh = sum(ashComposition)
    ashComposition = [speciesComp / totalAsh for speciesComp in ashComposition]

    moles = fuelMix.species_moles

    for index, species in enumerate(["C(gr)", "H", "O", "N", "S", "CL"]):
        speciesMass = ultimate[index] * (1 - ashFraction) * fuelDryMass
        speciesMoles = speciesMass / phases.Mw[species]
        moles[phases.indices[species]] = speciesMoles
    
    ashMass = ashFraction * fuelDryMass
    # To start, all ash is set to SiO2
    moles[phases.indices['SiO2(hqz)']] = ashMass / phases.Mw['SiO2(hqz)']
    fuelMix.species_moles = moles

    # Then, ash compostion is updated according to given distribution
    fuelMix.redistribute_ash(ashComposition)

    fuelMix.fuelHHV = HHV
    fuelMix.fuelLHV = LHV
    fuelMix.fuelAshFraction = ashFraction
    fuelMix.fuelDryMass = fuelDryMass
    fuelMix.set_moisture(moisture)

    return fuelMix

def create_air_stream(airMass, T=298.15, P=101325):
    airStream = phases.stream()
    airMoles = airMass / phases.Mw_air
    moles = airStream.species_moles
    moles[phases.indices['O2']] = 0.21*airMoles
    moles[phases.indices['N2']] = 0.78*airMoles
    moles[phases.indices['Ar']] = 0.01*airMoles
    airStream.species_moles = moles
    airStream.T = T
    airStream.P = P
    return airStream

def create_O2_stream(oxygenMass, T=298.15, P=101325):
    oxygenStream = phases.stream()
    moles = oxygenStream.species_moles
    moles[phases.indices['O2']] = oxygenMass / phases.Mw['O2']
    oxygenStream.species_moles = moles
    oxygenStream.T = T
    oxygenStream.P = P
    return oxygenStream

def create_steam_stream(steamMass, T=298.15, P=101325):
    steamStream = phases.stream()
    steamMoles = steamStream.species_moles
    steamMoles[phases.indices['H2O']] = steamMass / phases.Mw['H2O']
    steamStream.species_moles = steamMoles
    steamStream.T = T
    steamStream.P = P
    return steamStream

# # ------------------
# def determine_stoich_ash(fuelStream): # DEPRECATED???
#     totalMoles = sum(fuelStream.species_moles)

#     m = lambda e: fuelStream.element_moles(e)
#     stoicO2Mole = (m('C') + 0.25*m('H') - 0.5*m('O') + m('S') + 0.5*m('Cl') +
#                     + m('Si') + 0.5*m('Ca') + 0.75*m('Al') +
#                     + 0.75*m('Fe') + 0.25*m('Na') + 0.25*m('K') + 0.5*m('Mg') +
#                     + 1.25*m('P') + m('Ti') + 0.75*m('Cr'))
#     # EDITED 31-10-22: Changed some coefficients I thought were wrong (Nícolas)
#     stoichO2Ratio = stoicO2Mole / totalMoles

#     return stoichO2Ratio
# # ------------------

def determine_stoich(fuelStream):
    totalMoles = sum(fuelStream.species_moles)

    m = lambda sp: fuelStream.species_moles[phases.indices[sp]]
    stoicO2Mole = m('C(gr)') + 0.25*m('H') - 0.5*m('O') + m('S') + 0.5*m('CL')
    stoichO2Ratio = stoicO2Mole / totalMoles

    return stoichO2Ratio

def convert_air_ER(fuelStream, airMass):
    stoichO2Ratio = determine_stoich(fuelStream)
    stoichO2Moles = stoichO2Ratio * sum(fuelStream.species_moles)
    airMoles = airMass / phases.Mw_air
    currentO2Moles = airMoles * 0.21
    ER = currentO2Moles / stoichO2Moles
    return ER

def convert_ER_air(fuelStream, ER):
    stoichO2Ratio = determine_stoich(fuelStream)
    stoichO2Moles = stoichO2Ratio * sum(fuelStream.species_moles)
    currentO2Moles = stoichO2Moles * ER
    airMoles = currentO2Moles / 0.21
    airMass = airMoles * phases.Mw_air
    return airMass

def convert_O2_ER(fuelStream, oxygenMass):
    stoichO2Ratio = determine_stoich(fuelStream)
    stoichO2Moles = stoichO2Ratio * sum(fuelStream.species_moles)
    currentO2Moles = oxygenMass / phases.Mw['O2']
    ER = currentO2Moles / stoichO2Moles
    return ER

def convert_ER_O2(fuelStream, ER):
    stoichO2Ratio = determine_stoich(fuelStream)
    stoichO2Moles = stoichO2Ratio * sum(fuelStream.species_moles)
    currentO2Moles = stoichO2Moles * ER
    oxygenMass = currentO2Moles * phases.Mw['O2']
    return oxygenMass

def create_air_from_ER(fuelStream, ER):
    return create_air_stream(convert_ER_air(fuelStream, ER))

def create_O2_from_ER(fuelStream, ER):
    return create_O2_stream(convert_ER_O2(fuelStream, ER))





def combine_streams(stream1, stream2):
    stream1moles = stream1.species_moles
    stream2moles = stream2.species_moles
    combinedMoles = [stream1moles[i] + stream2moles[i] for i in range(len(stream1moles))]

    combinedStream = phases.stream()
    combinedStream.species_moles = combinedMoles

    return combinedStream

def blend_fuels(fuelStream1, fuelStream2):

    fuelBlend = phases.stream()

    return fuelBlend
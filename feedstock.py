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

def get_fuel_from_db(fuelID, fuelDryMass, moisture):
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

    fuelMix = phases.fuel()
    fuelMix.HHV = fu.get_hv(fuelID, 'HHV')
    fuelMix.LHV = fu.get_hv(fuelID, 'LHV')

    molesInOrder = np.zeros(len(fuelMix.species_names))
    for index, species in enumerate(fuelMix.species_names):
        if species in molesBySpecies.keys():
            molesInOrder[index] = molesBySpecies[species]
        else:
            molesInOrder[index] = 0

    fuelMix.species_moles = molesInOrder

    fuelMix.set_moisture(moisture)

    return fuelMix

def get_fuel_mix(fuelDryMass, ultimate, ashFraction, moisture, HHV, LHV=None, 
                ashComposition=[50, 20, 10, 10, 10, 0, 0, 0, 0, 0, 0]):
    '''
    Create a Cantera 'Mixture' object representing the given fuel.

    Parameters
    ----------
    fuelDryMass : float
        The fuel mass [kg] in dry basis
    ultimate : list
        Ultimate analysis of the fuel: [C H O N S Cl] (%m/m dry basis)
    ashFraction : float
        Fraction of ash in the fuel (%m/m dry basis)
    moisture : float
        Moisture content of the fuel (m/m fraction dry basis)
    HHV : float
        Higher heating value of the fuel (MJ/kg)
    LHV : float
        Lower heating value of the fuel (MJ/kg)

    Returns
    -------
    fuelMix : Fuel object
        The fuel mixture object.
    '''
    fuelMix = phases.fuel()

    # Normalization: ultimate + ashFraction should add up to 100%
    total = sum(ultimate) + ashFraction
    ultimate = [item / total for item in ultimate]
    ashFraction = ashFraction / total
    # Ash compositon should be normalized as well
    totalAsh = sum(ashComposition)
    ashComposition = [speciesComp / totalAsh for speciesComp in ashComposition]

    moles = fuelMix.species_moles

    for index, species in enumerate(["C(gr)", "H", "O", "N", "S", "CL"]):
        speciesMass = ultimate[index] * fuelDryMass
        speciesMoles = speciesMass / phases.Mw[species]
        moles[phases.indices[species]] = speciesMoles
    
    ashMass = ashFraction * fuelDryMass
    # To start, all ash is set to SiO2
    moles[phases.indices['SiO2(hqz)']] = ashMass / phases.Mw['SiO2(hqz)']
    fuelMix.species_moles = moles

    # Then, ash compostion is updated according to given distribution
    fuelMix.redistribute_ash(ashComposition)

    fuelMix.HHV = HHV
    fuelMix.LHV = LHV
    fuelMix.ashFraction = ashFraction
    fuelMix.set_moisture(moisture)

    return fuelMix

# TODO
def create_air_stream(airMass, T, P):
    airStream = phases.mix()
    return airStream

def create_O2_stream(O2Mass, T, P):
    oxygenStream = phases.mix()
    return oxygenStream

def create_steam_stream(steamMass, T, P):
    steamStream = phases.mix()
    return steamStream
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

def get_fuel_from_db(fuelID, fuelMass):
    '''
    Create a Cantera 'Mixture' object representing the given fuel.

    Parameters
    ----------
    fuelID : str
        The fuel ID as given by the CSV.
    fuelMass : float
        The fuel mass [kg]

    Returns
    -------
    fuelMix : Cantera 'Mixture' object
        The fuel mixture object.
    '''
    comp = fu.fuel_comp(fuelID)
    massBySpecies = {key: value*fuelMass for key, value in comp.items()} # kg

    molesBySpecies = {key: value/phases.get_molecular_weight(key) for key, value in massBySpecies.items()} # kmol

    fuelMix = phases.fuel()
    fuelMix.moisture = fu.moisture(fuelID)
    fuelMix.HHV = fu.get_hv(fuelID, 'HHV')
    fuelMix.LHV = fu.get_hv(fuelID, 'LHV')

    molesInOrder = np.zeros(len(fuelMix.species_names))
    for index, species in enumerate(fuelMix.species_names):
        if species in molesBySpecies.keys():
            molesInOrder[index] = molesBySpecies[species]
        else:
            molesInOrder[index] = 0

    fuelMix.species_moles = molesInOrder

    return fuelMix

def get_fuel_mix(ultimate, ashFraction, moisture, HHV, LHV=None):

    fuelMix = None
    return fuelMix
#!/usr/bin/env python

"""This script defines functions to get values from a fuel database.

The fuel database is a predefined CSV file ('fuels.csv'). The file can be 
edited using any spreadsheet editor to add new fuel compounds.

@author = Rodolfo Rodrigues
@contact = rodolfo.rodrigues@ufsm.br
@data = April, 2012, rev.: June, 2013 (adapted to use cython Cantera)
"""

# FIXME: 2013/05/25 Not possible to use a list of single value as parameter to 
# functions.
# FIXME: 2022/04/15 Ask professor if no data should return 0 or np.nan.
# TODO: 2022/04/27 Create functions: from fuelID and mass get Mix, from Mix get mass.
# TODO: 2022/04/18 Biochemical composition should be in csv.
# TODO: 2022/05/01 Find ash composition for petroleum coke.
# TODO: 2022/05/01 Add a cantera mixture for steam + water.
# TODO: 2022/05/01 Find heating values HHV and LHV from correlations.

# TODO LATER: 2022/05/01 Get heat of formation.

#==============================================================================
# import libraries
#==============================================================================
import numpy as np
import pandas as pd
import pp

# Create a dataframe from the csv
with open('fuels.csv','r') as f1:
    fuels = pd.read_csv(f1, sep=',', header=0, index_col=0)
    f1.close()

# create a list of the fuel names from the first column of the dataframe
myID = fuels.index.values

# Description, Type, Category, etc: read directly using fuels.loc
# Moisture, FC, VM: read directly then divide by 100
# HHV, LHV: read directly
# C, H, O, N, S, Cl: read directly then divide by 100

with open('ashes.csv', 'r') as f2:
    ashes = pd.read_csv(f2, sep=',', header=0, index_col=0)
    f2.close()

def ash(fuelID):
    '''
    Get the ash fraction value for a given fuel. 
    The fuel must be available in the database (file: 'fuels.csv').

    Parameters
    ----------
    fuelID : str
        The fuel ID as given by the CSV.

    Returns
    -------
    ashF : float
        Ash fraction [kg/kg]
    '''

    readAsh = fuels.loc[fuelID]['Ash']/100

    if pd.isnull(readAsh):
        ashF = 0
    else:
        ashF = readAsh

    return ashF


def ashComp(fuelID):
    '''
    Get the ash composition for a given fuel. 
    SiO2 CaO Al2O3 Fe2O3 Na2O K2O MgO P2O5 TiO2 SO3 Cr2O3
    The fuel must be available in the database (file: 'fuels.csv').
    
    Default values are used if they are not available. Those values are
    chosen from VASSILEV et al. (2013) according to type and category of fuel.

    Parameters
    ----------
    fuelID : str
        The fuel ID as given by the CSV.

    Returns
    -------
    comp : ndarray
        Composition of ash [kg/kg]
    '''
    ashF = ash(fuelID)
    rComp = fuels.loc[fuelID]['SiO2':'Cr2O3']/100

    compDF = 0

    # If ash fraction is 0, ash composition doesn't matter
    if ashF == 0:
        comp = np.ones(len(rComp))/len(rComp)

    else:
        comp = np.zeros(len(rComp))

    # If ash composition is not given by csv,
    # use values from Vassilev2013 given Type and Category
        if pd.isnull(rComp['SiO2']) and pd.isnull(rComp['CaO']):
            rType = fuels.loc[fuelID]['Type']
            rCat = fuels.loc[fuelID]['Category']
            if rType in ashes['Type'].values: 
                fType = rType
            else:
                fType = 'Other'
            if rCat in ashes['Category'].values: 
                fCat = rCat
            else:
                fCat = 'Other'
            compDF = ashes.loc[ashes['Type']==fType][ashes['Category']==fCat]
            comp = compDF.loc[:,'SiO2':'Cr2O3'].values[0]/100

    # If ash composition is given by csv, use it
        else:
            for index, species in enumerate(rComp):
                if pd.isnull(species):
                    comp[index] = 0
                else:
                    comp[index] = species

    return comp

def fuelComp(fuelID):
    '''
    Get the full mass composition for the fuel.
    The fuel must be available in the database (file: 'fuels.csv').

    Parameters
    ----------
    fuelID : str
        The fuel ID as given by the CSV.

    Returns
    -------
    fuelComp : ndarray
        An array representing the mass composition of each species.
    '''


def getFuelMix(fuelID, fuelMass):
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
    return 0